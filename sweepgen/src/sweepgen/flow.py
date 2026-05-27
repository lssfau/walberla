# This file is part of waLBerla. waLBerla is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# waLBerla is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.

from typing import Callable, overload, Iterable

from .sweep import Sweep
from pystencils import CreateKernelConfig
from pystencils.flow import FlowgraphNode
from pystencils.flow.builders import block, EquationsBlockBuilder

from pystencilssfg.composer import SfgBasicComposer


@overload
def sweep(
    *,
    config: CreateKernelConfig | None = None,
    preds: Iterable[FlowgraphNode] | None = None,
    name: str | None = None,
) -> Callable[[Callable[[EquationsBlockBuilder], None]], Sweep]: ...  # noqa: E704


@overload
def sweep(func: Callable[[EquationsBlockBuilder], None], /) -> Sweep: ...  # noqa: E704


def sweep(
    func: Callable[[EquationsBlockBuilder], None] | None = None,
    /,
    config: CreateKernelConfig | None = None,
    preds: Iterable[FlowgraphNode] | None = None,
    name: str | None = None,
):
    """Create a `Sweep` generator directly from a `pystencils.flow` equations block.

    This function behaves like `ps.flow.block <pystencils.flow.block>`,
    but in addition immediately creates a `Sweep` code generator object from the given equations block.

    **Example:** The following example creates the generator for a sweep called ``FieldCopy``
    copying all values from one field ``g`` to another field ``f``,
    and then invokes code generation via ``sfg.generate``:

    .. code-block:: Python

        import sweepgen as sg

        ...

        cfg = ps.CreateKernelConfig()

        @sg.flow.sweep(config=cfg)
        def FieldCopy(_eq):
            _eq.store[f()] = g()

        sfg.generate(FieldCopy)

    Args:
        config: Code generation config for this sweep
        preds: Predecessors to this ``ps.flow`` block; see `ps.flow.block <pystencils.flow.block>`
        name: Name for this sweep; if none is given, the decorated function's name is used
    """

    def decorate(func: Callable[[EquationsBlockBuilder], None]) -> Sweep:
        eqblock = block(preds=preds, name=name)(func)
        return Sweep(eqblock.name, eqblock, config=config)

    if func is None:
        return decorate
    else:
        return decorate(func)


def generate_sweep(
    sfg: SfgBasicComposer,
    *,
    config: CreateKernelConfig | None = None,
    preds: Iterable[FlowgraphNode] | None = None,
    name: str | None = None,
) -> Callable[[Callable[[EquationsBlockBuilder], None]], None]:
    """Generate a sweep directly from a `pystencils.flow` equations block.

    This function behaves like `ps.flow.block <pystencils.flow.block>`,
    but in addition immediately generates code for a sweep from the given equations block.

    **Example:** The following example generates a sweep called ``FieldCopy`` copying
    all values from one field ``g`` to another field ``f``:

    .. code-block:: Python

        import sweepgen as sg

        ...

        cfg = ps.CreateKernelConfig()

        @sg.flow.generate_sweep(sfg, config=cfg)
        def FieldCopy(_eq):
            _eq.store[f()] = g()

    .. admonition:: Example Apps
        :class: seealso

        - :walberla-example-app:`DoubleShearLayer`
        - :walberla-example-app:`ParallelPlates`

    Args:
        sfg: Active SFG composer object
        config: Code generation config for this sweep
        preds: Predecessors to this ``ps.flow`` block; see `ps.flow.block <pystencils.flow.block>`
        name: Name for this sweep; if none is given, the decorated function's name is used
    """

    def decorate(func: Callable[[EquationsBlockBuilder], None]) -> None:
        sw = sweep(config=config, preds=preds, name=name)(func)
        sfg.generate(sw)

    return decorate
