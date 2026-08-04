# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2025 Frederik Hennig <frederik.hennig@fau.de>

from __future__ import annotations

from typing import Sequence, Iterable
from dataclasses import dataclass
from itertools import chain

from pystencilssfg.composer.class_composer import SfgMethodSequencer

from pystencils.types import PsType, constify, deconstify

from pystencilssfg import SfgComposer
from pystencilssfg.lang import (
    VarLike,
    asvar,
    SfgVar,
    AugExpr,
    Ref,
    CppClass,
    cpptype,
    ExprLike,
)
from pystencilssfg.ir import SfgSequence

from ..api import (
    StructuredBlockForest,
    SharedPtr,
    SweepConstraintError,
)


@dataclass(frozen=True)
class PropertyConstraint:
    """Constraint on a property's values"""

    condition: AugExpr
    message: str

    @staticmethod
    def equals(a: ExprLike, b: ExprLike, msg: str):
        return PropertyConstraint(AugExpr.format("{} == {}", a, b), msg)

    def render(self, sfg: SfgComposer) -> SfgSequence:
        return sfg.seq(
            sfg.branch(sfg.expr("!({})", self.condition))(
                sfg.expr("throw {};", SweepConstraintError().ctor(f'"{self.message}"'))
            )
        )


@dataclass(frozen=True)
class Property:
    name: str
    dtype: PsType
    getter: bool
    setter: bool
    initializer: tuple[AugExpr, ...] | None = None
    constraints: tuple[PropertyConstraint, ...] = ()

    def __post_init__(self):
        if self.dtype.const and self.setter:
            raise ValueError("Constant properties cannot have setters.")

    @property
    def var(self) -> SfgVar:
        return SfgVar(self.name, self.dtype)

    @property
    def member_var(self) -> SfgVar:
        return SfgVar(self.name + "_", self.dtype)


class PropertiesContainer(CppClass):
    properties: dict[str, Property]
    ctor_params: tuple[SfgVar, ...]

    def ctor(self, **kwargs):
        if set(kwargs.keys()) != set(p.name for p in self.ctor_params):
            raise ValueError(f"Invalid parameter list to constructor of {self.dtype}")

        return self.ctor_bind(*(kwargs[p.name] for p in self.ctor_params))

    def get(self, prop: str | VarLike) -> AugExpr:
        if not isinstance(prop, str):
            prop = asvar(prop)
            prop_name = prop.name
        else:
            prop_name = prop

        if prop_name not in self.properties:
            raise ValueError(f"Invalid property: {prop_name}")

        return AugExpr.format("{}.{}_", self, prop_name)

    def render_forwarding_ctor(
        self, sfg: SfgComposer
    ) -> SfgComposer.ConstructorBuilder:
        ctor = (
            sfg.constructor(*self.ctor_params)
            .init(self)(*self.ctor_params)
            .body(
                *(
                    chain.from_iterable(
                        (
                            (constr.render(sfg) for constr in prop.constraints)
                            for prop in self.properties.values()
                        )
                    )
                )
            )
        )
        return ctor

    def render_public_interface(
        self, sfg: SfgComposer
    ) -> tuple[SfgMethodSequencer, ...]:
        def make_methods(p: Property) -> Sequence[SfgMethodSequencer]:
            methods = []

            if p.getter:
                methods.append(
                    sfg.method(f"{p.name}")
                    .returns(Ref(constify(p.dtype)))
                    .const()
                    .inline()(f"return {self}.{p.name}_;")
                )

            if p.setter:
                methods.append(
                    sfg.method(f"{p.name}")
                    .returns(Ref(p.dtype))
                    .inline()(f"return {self}.{p.name}_;")
                )

            return methods

        methods = chain.from_iterable(make_methods(p) for p in self.properties.values())
        return tuple(methods)


class PropertiesContainerBuilder:
    """Collect and manage properties of sweep classes."""

    blockforest_shared_ptr = StructuredBlockForest.shared_ptr()

    def __init__(self) -> None:
        self._properties: dict[str, Property] = dict()

    @property
    def properties(self) -> Iterable[Property]:
        return self._properties.values()

    def add_property(
        self,
        prop: VarLike,
        const=False,
        setter: bool = True,
        getter: bool = True,
        initializer: tuple[AugExpr, ...] | None = None,
        constraints: tuple[PropertyConstraint, ...] = (),
    ):
        prop_var = asvar(prop)
        if prop_var.name in self._properties:
            raise ValueError(f"Duplicate property: {prop_var}")

        dtype = constify(prop_var.dtype) if const else deconstify(prop_var.dtype)

        self._properties[prop_var.name] = Property(
            prop_var.name, dtype, getter, setter, initializer, constraints
        )

    def add_blockforest_shared_ptr(self) -> SharedPtr:
        self.add_property(
            self.blockforest_shared_ptr,
            getter=False,
            setter=False,
            initializer=(StructuredBlockForest.shared_ptr_ref(),),
        )
        return self.blockforest_shared_ptr

    def render_struct(self, sfg: SfgComposer, name: str) -> type[PropertiesContainer]:
        ctor = sfg.constructor()

        for p in self.properties:
            if p.initializer is None:
                ctor.add_param(p.var)
            else:
                for var in chain.from_iterable(e.depends for e in p.initializer):
                    if var not in ctor.parameters:
                        if var.name == asvar(self.blockforest_shared_ptr).name:
                            ctor.add_param(var, 0)
                        else:
                            ctor.add_param(var)

        for p in self._properties.values():
            if p.initializer is not None:
                ctor.init(p.member_var)(*p.initializer)
            else:
                ctor.init(p.member_var)(p.var)

        sfg.struct(name)(
            #   Members
            *(p.member_var for p in self.properties),
            #   Constructor
            ctor,
        )

        namespace = sfg._cursor.current_namespace.fqname

        class Container(PropertiesContainer):
            properties = self._properties.copy()
            ctor_params = tuple(ctor.parameters)
            template = cpptype(f"{namespace}::{name}")

        return Container
