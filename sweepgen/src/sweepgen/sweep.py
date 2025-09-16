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

from typing import Sequence, cast
from dataclasses import dataclass
from collections import defaultdict
from abc import ABC, abstractmethod
import warnings

from pystencilssfg.composer.custom import CustomGenerator

from pystencils import (
    Assignment,
    AssignmentCollection,
    CreateKernelConfig,
    Field,
    FieldType,
    Target,
    create_type,
)
from pystencils.types import deconstify, PsCustomType, PsStructType

from pystencilssfg import SfgComposer
from pystencilssfg.lang import (
    asvar,
    SfgVar,
    AugExpr,
    SupportsFieldExtraction,
    SupportsVectorExtraction,
)
from pystencilssfg.lang.cpp import std
from pystencilssfg.ir import SfgKernelHandle, SfgCallTreeNode
from .build_config import WalberlaBuildConfig
from .core.reflection import GeneratedClassWrapperBase
from .api import (
    StructuredBlockForest,
    GenericWalberlaField,
    GhostLayerFieldPtr,
    GpuFieldPtr,
    IndexListBufferPtr,
    SparseIndexList,
    IBlockPtr,
    BlockDataID,
    Vector2,
    Vector3,
    CellInterval,
    CellIdx,
)

from .core.properties import PropertiesContainerBuilder, PropertiesContainer
from .core.blockforest_extraction import BlockforestParamExtraction
from .core.suppress_diagnostics import SuppressDiagnostics


@dataclass
class FieldInfo(ABC):
    field: Field

    @property
    @abstractmethod
    def entity(self) -> AugExpr: ...

    @property
    @abstractmethod
    def view(self) -> AugExpr: ...

    @abstractmethod
    def create_view(self, block: IBlockPtr) -> AugExpr: ...

    def view_extraction(self) -> SupportsFieldExtraction:
        assert isinstance(self.view, SupportsFieldExtraction)
        return self.view


@dataclass
class GlFieldInfo(FieldInfo):
    glfield: GenericWalberlaField
    data_id: BlockDataID

    @property
    def entity(self) -> AugExpr:
        return self.data_id

    @property
    def view(self) -> AugExpr:
        return self.glfield

    def create_view(self, block: IBlockPtr) -> AugExpr:
        return block.getData(self.glfield.field_type, self.data_id)


@dataclass
class IndexListInfo(FieldInfo):
    idx_list: SparseIndexList
    idx_vector_view: std.span

    @property
    def entity(self) -> AugExpr:
        return self.idx_list

    @property
    def view(self) -> AugExpr:
        return self.idx_vector_view

    def create_view(self, block: IBlockPtr) -> AugExpr:
        return self.idx_list.view(block.deref())


class ShadowFieldsManager:
    """
    Information needed:
        - Name and data type of field
        - Name of shadow field
    """

    def __init__(self) -> None:
        self._originals: dict[str, GlFieldInfo] = dict()

    def add_shadow(self, original: GlFieldInfo):
        self._originals[original.field.name] = original

    def get_shadow(self, original_name: str) -> AugExpr:
        original = self._originals[original_name]
        return AugExpr.format(
            "this->{}({})", self._getter(original.field.name), original.glfield
        )

    def perform_swap(self, original_name: str, shadow_field: GlFieldInfo) -> AugExpr:
        original = self._originals[original_name]
        return original.glfield.swapDataPointers(shadow_field.glfield)

    def render(self, sfg: SfgComposer):
        if not self._originals:
            return ()

        sfg.include("<memory>")

        cache_ptrs = []
        getters = []

        for orig_name, orig in self._originals.items():
            unique_ptr_type = PsCustomType(
                f"std::unique_ptr< {orig.glfield.field_type} >"
            )
            cache_ptr_name = self._cache_ptr(orig.glfield)
            cache_ptrs.append(sfg.var(cache_ptr_name, unique_ptr_type))

            getters.append(
                sfg.method(self._getter(orig_name))
                .returns(orig.glfield.get_dtype())
                .inline()(
                    sfg.branch(f"{cache_ptr_name} == nullptr")(
                        AugExpr.format(
                            "{}.reset({});",
                            cache_ptr_name,
                            orig.glfield.cloneUninitialized(),
                        )
                    ),
                    f"return {cache_ptr_name}.get();",
                )
            )

        return (sfg.private(*cache_ptrs, *getters),)

    def _getter(self, orig_name: str) -> str:
        return f"_getShadow_{orig_name}"

    def _cache_ptr(self, orig: GenericWalberlaField) -> str:
        return f"shadow_{str(orig)}_"


def combine_vectors(
    scalars: set[SfgVar],
) -> dict[SupportsVectorExtraction, tuple[SfgVar, ...]]:
    """Attempt to combine vector component symbols into vectors.

    This function modifies the `scalars` parameter in-place by removing grouped scalar
    components and replacing them by their vectors.
    """
    potential_vectors: defaultdict[str, list[tuple[int, SfgVar]]] = defaultdict(list)

    for s in scalars:
        tokens = s.name.split("_")
        try:
            coord = int(tokens[-1])
            vec_name = "_".join(tokens[:-1])
            potential_vectors[vec_name].append((coord, s))
        except ValueError:
            pass

    all_vectors: dict[SupportsVectorExtraction, tuple[SfgVar, ...]] = dict()

    for name, entries in potential_vectors.items():
        entries = sorted(entries, key=lambda t: t[0])
        coords = [e[0] for e in entries]
        components = tuple(e[1] for e in entries)
        scalar_types = set(e[1].dtype for e in entries)

        vec_type: type[AugExpr]
        if set(coords) == {0, 1}:
            vec_type = Vector2
        elif set(coords) == {0, 1, 2}:
            vec_type = Vector3
        else:
            continue

        if len(scalar_types) > 1:
            warnings.warn(
                f"Skipping vector grouping for vector {name}: "
                f"Encountered vector components {components} with conflicting data types {scalar_types}.",
                UserWarning,
            )
            continue

        name_conflict = False
        for s in scalars:
            if s.name == name:
                warnings.warn(
                    f"Skipping vector grouping for vector {name}: "
                    f"Vector name conflicts with existing scalar parameter {s}.",
                    UserWarning,
                )
                name_conflict = True

        if name_conflict:
            continue

        scalar_type = deconstify(scalar_types.pop())
        vector = vec_type(scalar_type).var(name)
        all_vectors[vector] = components

        for c in components:
            scalars.remove(c)
        scalars.add(asvar(vector))

    return all_vectors


class Sweep(CustomGenerator):
    """Generate a waLBerla sweep from a pystencils kernel.

    Args:
        name: Name of the sweep class
        assignments: Assignment collection defining the pystencils kernel
        config: Optional code generator configuration. Options specified in this configuration
            override options from the `WalberlaBuildConfig`.
    """

    def __init__(
        self,
        name: str,
        assignments: Sequence[Assignment] | AssignmentCollection,
        config: CreateKernelConfig | None = None,
    ):
        if config is not None:
            cfg = config.copy()

            if cfg.get_option("ghost_layers") is not None:
                raise ValueError(
                    "Specifying `ghost_layers` in your codegen config is invalid when generating a waLBerla sweep."
                )

            if (
                cfg.get_option("iteration_slice") is None
                and cfg.get_option("index_field") is None
            ):
                cfg.ghost_layers = 0
        else:
            cfg = CreateKernelConfig(ghost_layers=0)

        manual_grid: bool = cfg.gpu.get_option("manual_launch_grid")
        if manual_grid:
            raise ValueError(
                "Setting `gpu.manual_launch_grid = True` is invalid for waLBerla sweeps."
            )

        self._name = name
        self._gen_config = cfg

        if isinstance(assignments, AssignmentCollection):
            self._assignments = assignments
        else:
            self._assignments = AssignmentCollection(assignments)  # type: ignore

        #   Set only later once the full codegen config is known
        self._glfield_type: type[GpuFieldPtr] | type[GhostLayerFieldPtr]

        #   Map from shadow field to shadowed field
        self._shadow_fields: dict[Field, Field] = dict()

        #   RESULTS - unset at this point
        self._generated_class: type[GeneratedClassWrapperBase] | None = None
        self._properties_container: type[PropertiesContainer] | None = None

    #   CONFIGURATION

    @property
    def sparse(self) -> bool:
        """If set to ``True``, the sweep is generated as a sparse sweep.

        Instead of iterating over all cells on a block,
        sparse sweeps process only cells from a given *index vector*.
        They can be used to implement phenomena that should only be applied
        in small regions of the simulation space,
        such as boundary conditions.

        If the ``config`` argument passed to the constructor has an
        `index_field <pystencils.codegen.config.CreateKernelConfig.index_field>` set,
        the sweep is automatically marked *sparse*.
        """
        return self._gen_config.get_option("index_field") is not None

    @sparse.setter
    def sparse(self, sparse_iteration: bool):
        if sparse_iteration:
            if self._gen_config.get_option("index_field") is not None:
                return

            if self._gen_config.get_option("index_dtype") != create_type("int64"):
                raise ValueError(
                    "Sparse sweeps require `int64_t` as their index data type. Check your code generator config."
                )

            self._gen_config.index_field = Field.create_generic(
                "indexList", 1, CellIdx, field_type=FieldType.INDEXED
            )
            self._gen_config.ghost_layers = None
        else:
            self._gen_config.index_field = None
            self._gen_config.ghost_layers = 0

    def swap_fields(self, field: Field, shadow_field: Field):
        """Register two fields as a swapping pair.

        The generated sweep will treat ``shadow_field`` as a temporary shadow of ``field``.
        The shadow field will not occur in the sweep's public API,
        but me managed internally.
        The data arrays of ``field`` and ``shadow_field`` will be swapped
        after each invocation of the sweep, such that data written to the shadow field
        overwrites the old content of ``field``.
        """
        if field in self._shadow_fields:
            raise ValueError(f"Field swap for {field} was already registered.")
        if shadow_field in self._shadow_fields:
            raise ValueError(f"Field swap for {shadow_field} was already registered.")
        if field.dtype != shadow_field.dtype:
            raise ValueError(
                "Field and its shadow must have the same element type for swapping"
            )

        self._shadow_fields[shadow_field] = field

    #   CODE GENERATION

    def _set_field_interface(self, target: Target):
        if target.is_gpu():
            self._glfield_type = GpuFieldPtr
        elif target.is_cpu():
            self._glfield_type = GhostLayerFieldPtr
        else:
            raise ValueError(
                f"Cannot generate sweep for target {self._gen_config.target}"
            )

    def _walberla_field(self, f: Field) -> GenericWalberlaField | IndexListBufferPtr:
        match f.field_type:
            case FieldType.GENERIC | FieldType.CUSTOM:
                return self._glfield_type.create(f)
            case FieldType.INDEXED:
                assert isinstance(f.dtype, PsStructType)
                return IndexListBufferPtr(f.dtype).var(f.name)
            case _:
                raise ValueError(
                    f"Unable to map field {f} of type {f.field_type} to a waLBerla field."
                )

    def _make_field_info(self, f: Field, target: Target) -> FieldInfo:
        match f.field_type:
            case FieldType.GENERIC | FieldType.CUSTOM:
                glfield = self._glfield_type.create(f)
                data_id = BlockDataID().var(f"{f.name}Id")
                return GlFieldInfo(f, glfield, data_id)
            case FieldType.INDEXED:
                assert isinstance(f.dtype, PsStructType)
                idx_list = SparseIndexList.from_field(f, target=target)
                view = idx_list.view_type().var(f"{f.name}_view")
                return IndexListInfo(f, idx_list, view)
            case _:
                raise ValueError(f"Unexpected field type: {f.field_type} at field  {f}")

    def _render_invocation(
        self, sfg: SfgComposer, target: Target, khandle: SfgKernelHandle
    ) -> tuple[SfgCallTreeNode, set[SfgVar]]:
        """Render and return the kernel invocation plus a set of additional parameters required
        at the call site."""

        if target.is_gpu():
            # from pystencils.codegen.config import GpuIndexingScheme

            #   TODO: Want default values for properties first,
            #   to define default stream and block size values
            # indexing_scheme = self._gen_config.gpu.get_option("indexing_scheme")
            # if indexing_scheme == GpuIndexingScheme.Linear3D:
            #     block_size = sfg.gpu_api.dim3(const=True).var("gpuBlockSize")
            #     return (sfg.gpu_invoke(khandle, block_size=block_size), {block_size})
            # else:
            return (sfg.gpu_invoke(khandle), set())

        else:
            return (sfg.call(khandle), set())

    def generate(self, sfg: SfgComposer) -> None:
        build_config = WalberlaBuildConfig.from_sfg(sfg)
        gen_config = build_config.get_pystencils_config()
        gen_config.override(self._gen_config)

        target = gen_config.get_target()
        self._set_field_interface(target)

        assignments = BlockforestParamExtraction.process(self._assignments)

        with SuppressDiagnostics(sfg, "unused-variable"):
            knamespace = sfg.kernel_namespace(f"{self._name}_kernels")
            khandle = knamespace.create(assignments, self._name, gen_config)

        ker_invocation, ker_call_site_params = self._render_invocation(
            sfg, target, khandle
        )

        all_fields: dict[str, FieldInfo] = {
            f.name: self._make_field_info(f, target) for f in khandle.fields
        }

        swaps: dict[str, GlFieldInfo] = dict()
        shadows_cache = ShadowFieldsManager()

        for shadow, orig in self._shadow_fields.items():
            shadow_fi = all_fields[shadow.name]
            assert isinstance(shadow_fi, GlFieldInfo)
            original_fi = all_fields[orig.name]
            assert isinstance(original_fi, GlFieldInfo)
            swaps[orig.name] = shadow_fi
            shadows_cache.add_shadow(original_fi)

        block = IBlockPtr().var("block")
        block_fields = sorted(
            (fi for fi in all_fields.values() if fi.field not in self._shadow_fields),
            key=lambda fi: fi.field.name,
        )

        props_builder = PropertiesContainerBuilder()

        parameters = khandle.scalar_parameters | ker_call_site_params

        blockforest_ref = StructuredBlockForest(ref=True).var("blocks")
        blockforest_params = BlockforestParamExtraction(blockforest_ref, block)
        parameters, need_blockforest = blockforest_params.filter_params(parameters)

        if need_blockforest:
            blockforest_ptr = props_builder.add_blockforest_shared_ptr()

        vector_groups = combine_vectors(parameters)

        for fi in block_fields:
            props_builder.add_property(fi.entity, setter=False, getter=True)

        for s in sorted(parameters, key=lambda p: p.name):
            props_builder.add_property(s, setter=True, getter=True)

        with sfg.namespace("detail"):
            props_struct = props_builder.render_struct(sfg, self._name + "Properties")

        property_cache = props_struct().var("properties_")
        property_cache_ref = props_struct(ref=True).bind(f"this->{property_cache}")

        def render_runmethod(ci: CellInterval | None = None):
            return [
                #   Get IDs from class
                *(
                    sfg.init(fi.entity)(property_cache_ref.get(fi.entity))
                    for fi in block_fields
                ),
                #   Extract field views
                *(sfg.init(fi.view)(fi.create_view(block)) for fi in block_fields),
                #   Get shadow fields from cache
                *(
                    sfg.init(shadow_info.glfield)(shadows_cache.get_shadow(orig_name))
                    for orig_name, shadow_info in swaps.items()
                ),
                #   Map GhostLayerFields to pystencils fields
                *(
                    sfg.map_field(fi.field, fi.glfield.with_cell_interval(ci))
                    for fi in all_fields.values()
                    if isinstance(fi, GlFieldInfo)
                ),
                #   Extract indexing information from field view for all remaining fields
                *(
                    sfg.map_field(fi.field, fi.view_extraction())
                    for fi in all_fields.values()
                    if not isinstance(fi, GlFieldInfo)
                ),
                #   Get parameters from class
                *(
                    sfg.init(param)(property_cache_ref.get(param))
                    for param in parameters
                ),
                #   Extract scalars from vectors
                *(
                    sfg.map_vector(components, vector)
                    for vector, components in vector_groups.items()
                ),
                #   Extract geometry information
                *(
                    (
                        sfg.init(blockforest_ref)(
                            AugExpr.format(
                                "*({})", property_cache_ref.get(blockforest_ptr)
                            )
                        ),
                    )
                    if need_blockforest
                    else ()
                ),
                *(blockforest_params.render_extractions(sfg, ci)),
                #   Invoke the kernel
                ker_invocation,
                #   Perform field swaps
                *(
                    shadows_cache.perform_swap(orig_name, shadow_info)
                    for orig_name, shadow_info in swaps.items()
                ),
            ]

        runmethods = [sfg.method("operator()")(*render_runmethod())]

        if not self.sparse:
            runmethods.append(
                sfg.method("runOnCellInterval")(
                    *render_runmethod(CellInterval(const=True, ref=True).var("ci"))
                )
            )

        sfg.klass(self._name)(
            sfg.public(
                property_cache.render_forwarding_ctor(sfg),
                *property_cache.render_public_interface(sfg),
                *runmethods,
            ),
            sfg.private(property_cache),
            *shadows_cache.render(sfg),
        )

        gen_class = sfg._cursor.get_entity(self._name)
        assert gen_class is not None
        from pystencilssfg.ir.entities import SfgClass

        assert isinstance(gen_class, SfgClass)

        class GenClassWrapper(GeneratedClassWrapperBase):
            _class = cast(SfgClass, gen_class)

        self._generated_class = GenClassWrapper
        self._properties_container = props_struct

    #   CODE-GENERATION RESULTS
    #   These properties are only available after the sweep was generated

    @property
    def generated_class(self) -> type[GeneratedClassWrapperBase]:
        if self._generated_class is None:
            raise AttributeError(
                "Generated class is unavailable - code generation was not run yet."
            )

        return self._generated_class

    @property
    def properties_container(self) -> type[PropertiesContainer]:
        if self._properties_container is None:
            raise AttributeError(
                "Properties container is unavailable - code generation was not run yet."
            )

        return self._properties_container
