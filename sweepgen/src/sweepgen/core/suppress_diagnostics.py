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


from __future__ import annotations

from enum import Enum, auto

from pystencilssfg import SfgComposer


class Diagnostic(Enum):
    UNUSED_VARIABLE = auto()
    CONVERSION = auto()


class SuppressDiagnostics:
    """Context manager that suppresses compiler diagnostics.

    .. warning::

        HIGHLY EXPERIMENTAL - GCC ONLY!
        This class is an experimental prototype and shall be extended in the near future
        to support further diagnostics and a wider range of compilers.
        See also walberla/walberla#250.
    """

    #   See https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html

    def __init__(self, sfg: SfgComposer, *diags: str | Diagnostic, impl: bool = True):
        self._sfg = sfg

        self._diagnostics: set[Diagnostic] = set()
        self._impl = impl

        self.suppress(*diags)

    def suppress(self, *diags: str | Diagnostic) -> SuppressDiagnostics:
        for diag in diags:
            if isinstance(diag, str):
                diag = Diagnostic[diag.upper().replace("-", "_")]
            self._diagnostics.add(diag)
        return self

    def __enter__(self):
        if self._diagnostics:
            self._sfg.code("#pragma GCC diagnostic push", impl=self._impl)

            for diag in self._diagnostics:
                diag_text = "-W" + diag.name.lower().replace("_", "-")
                self._sfg.code(f'#pragma GCC diagnostic ignored "{diag_text}"', impl=self._impl)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self._sfg.code("#pragma GCC diagnostic pop", impl=self._impl)
