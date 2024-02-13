from pystencils import Target
from pystencils.stencil import inverse_direction
from pystencils.typing import BasicType

from lbmpy.advanced_streaming import AccessPdfValues, numeric_offsets, numeric_index, Timestep, is_inplace
# until lbmpy version 1.3.2 
try:
    from lbmpy.advanced_streaming.indexing import MirroredStencilDirections
except ImportError:
    from lbmpy.custom_code_nodes import MirroredStencilDirections
from lbmpy.boundaries.boundaryconditions import LbBoundary
from lbmpy.boundaries import (ExtrapolationOutflow, FreeSlip, UBB, DiffusionDirichlet,
                              NoSlipLinearBouzidi, QuadraticBounceBack)

from pystencils_walberla.additional_data_handler import AdditionalDataHandler


interpolation_bc_check_template = """
if(!isFlagSet(it.neighbor({cx}, {cy}, {cz}, 0), domainFlag)){{ 
   //Linear-Bouzidi requires 2 fluid nodes: if the 2nd node is not available abort, 
   //apply Bounce Back at that point. This clearly lowers the accuracy and makes inconsistent the
   //calculation of the total force
   element.q = -1.0;
   WALBERLA_LOG_INFO_ON_ROOT("Warning: Bouzidi cannot be applied at least on one boundary link.")
}} //end if to check Bouzidi applicability  
"""


def default_additional_data_handler(boundary_obj: LbBoundary, lb_method, field_name, target=Target.CPU,
                                    pdfs_data_type=None, zeroth_timestep=None):
    if not boundary_obj.additional_data:
        return None
    if isinstance(boundary_obj, FreeSlip):
        return FreeSlipAdditionalDataHandler(lb_method.stencil, boundary_obj)
    elif isinstance(boundary_obj, UBB):
        return UBBAdditionalDataHandler(lb_method.stencil, boundary_obj)
    elif isinstance(boundary_obj, ExtrapolationOutflow):
        return OutflowAdditionalDataHandler(lb_method.stencil, boundary_obj, target=target, field_name=field_name,
                                            pdfs_data_type=pdfs_data_type, zeroth_timestep=zeroth_timestep)
    elif isinstance(boundary_obj, NoSlipLinearBouzidi):
        return NoSlipLinearBouzidiAdditionalDataHandler(lb_method.stencil, boundary_obj)
    elif isinstance(boundary_obj, QuadraticBounceBack):
        return QuadraticBounceBackAdditionalDataHandler(lb_method.stencil, boundary_obj)
    else:
        raise ValueError(f"No default AdditionalDataHandler available for boundary of type {boundary_obj.__class__}")


class FreeSlipAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, FreeSlip)
        super(FreeSlipAdditionalDataHandler, self).__init__(stencil=stencil)

    def data_initialisation(self, direction):
        def array_pattern(dtype, name, content):
            return f"const {str(dtype)} {name} [] = {{ {','.join(str(c) for c in content)} }};"

        offset = self._walberla_stencil[direction]
        inv_offset = inverse_direction(self._walberla_stencil[direction])
        offset_dtype = "int32_t"
        mirror_stencil = MirroredStencilDirections.mirror_stencil

        axis_mirrored = []
        for i, name in enumerate(["x", "y", "z"]):
            mirrored_dir = [self._walberla_stencil.index(mirror_stencil(d, i)) for d in self._walberla_stencil]
            axis_mirrored.append(array_pattern(offset_dtype, f"{name}_axis_mirrored_stencil_dir", mirrored_dir))

        init_list = axis_mirrored[0:self._dim]

        init_list += [f"const Cell n = it.cell() + Cell({offset[0]}, {offset[1]}, {offset[2]});",
                      f"int32_t ref_dir = {self._walberla_stencil.index(inv_offset)}; // dir: {direction}",
                      "element.wnx = 0; // compute discrete normal vector of free slip wall",
                      "element.wny = 0;",
                      f"if( flagField->isPartOfMaskSet( n.x() + {inv_offset[0]}, n.y(), n.z(), domainFlag ) )",
                      "{",
                      f"   element.wnx = {inv_offset[0]};",
                      "   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];",
                      "}",
                      f"if( flagField->isPartOfMaskSet( n.x(), n.y() + {inv_offset[1]}, n.z(), domainFlag ) )",
                      "{",
                      f"   element.wny = {inv_offset[1]};",
                      "   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];",
                      "}"]
        if self._dim == 3:
            init_list += ["element.wnz = 0;",
                          f"if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + {inv_offset[2]}, domainFlag ) )",
                          "{",
                          f"   element.wnz = {inv_offset[2]};",
                          "   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];",
                          "}",
                          "// concave corner (neighbors are non-fluid)",
                          "if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )",
                          "{",
                          f"   element.wnx = {inv_offset[0]};",
                          f"   element.wny = {inv_offset[1]};",
                          f"   element.wnz = {inv_offset[2]};",
                          f"   ref_dir = {self._walberla_stencil.index(inv_offset)};",
                          "}"]
        elif self._dim == 2:
            init_list += ["// concave corner (neighbors are non-fluid)",
                          "if( element.wnx == 0 && element.wny == 0 )",
                          "{",
                          f"   element.wnx = {inv_offset[0]};",
                          f"   element.wny = {inv_offset[1]};",
                          f"   ref_dir = {self._walberla_stencil.index(inv_offset)};",
                          "}"]
        init_list.append("element.ref_dir = ref_dir;")

        return "\n".join(init_list)


class UBBAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, UBB)
        super(UBBAdditionalDataHandler, self).__init__(stencil=stencil)

    @property
    def constructor_argument_name(self):
        return "velocityCallback"

    @property
    def constructor_arguments(self):
        return f", std::function<Vector3<real_t>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& " \
               f"{self.constructor_argument_name} "

    @property
    def initialiser_list(self):
        return "elementInitialiser(velocityCallback),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, *_):
        init_list = ["Vector3<real_t> InitialisationAdditionalData = elementInitialiser(Cell(it.x(), it.y(), it.z()), "
                     "blocks, *block);", "element.vel_0 = InitialisationAdditionalData[0];",
                     "element.vel_1 = InitialisationAdditionalData[1];"]
        if self._dim == 3:
            init_list.append("element.vel_2 = InitialisationAdditionalData[2];")

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return "std::function<Vector3<real_t>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> " \
               "elementInitialiser; "


class NoSlipLinearBouzidiAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, NoSlipLinearBouzidi)

        self._dtype = BasicType(boundary_object.data_type).c_name
        self._blocks = "const shared_ptr<StructuredBlockForest>&, IBlock&)>"
        super(NoSlipLinearBouzidiAdditionalDataHandler, self).__init__(stencil=stencil)

    @property
    def constructor_argument_name(self):
        return "wallDistanceBouzidi"

    @property
    def constructor_arguments(self):
        return f", std::function<{self._dtype}(const Cell &, const Cell &, {self._blocks}&" \
               f"{self.constructor_argument_name} "

    @property
    def initialiser_list(self):
        return f"elementInitialiser({self.constructor_argument_name}),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, direction):
        cx = self._walberla_stencil[direction][0]
        cy = self._walberla_stencil[direction][1]
        cz = self._walberla_stencil[direction][2]
        fluid_cell = "Cell(it.x(), it.y(), it.z())"
        boundary_cell = f"Cell(it.x() + {cx}, it.y() + {cy}, it.z() + {cz})"
        check_str = interpolation_bc_check_template.format(cx=-cx, cy=-cy, cz=-cz)
        init_element = f"elementInitialiser({fluid_cell}, {boundary_cell}, blocks, *block)"
        init_list = [f"const {self._dtype} q = (({self._dtype}) {init_element});",
                     "element.q = q;",
                     check_str]

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return f"std::function<{self._dtype}(const Cell &, const Cell &, {self._blocks} elementInitialiser; "


class QuadraticBounceBackAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, QuadraticBounceBack)

        self._dtype = BasicType(boundary_object.data_type).c_name
        self._blocks = "const shared_ptr<StructuredBlockForest>&, IBlock&)>"
        super(QuadraticBounceBackAdditionalDataHandler, self).__init__(stencil=stencil)

    @property
    def constructor_argument_name(self):
        return "wallDistanceQuadraticBB"

    @property
    def constructor_arguments(self):
        return f", std::function<{self._dtype}(const Cell &, const Cell &, {self._blocks}&" \
               f"{self.constructor_argument_name} "

    @property
    def initialiser_list(self):
        return f"elementInitialiser({self.constructor_argument_name}),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, direction):
        cx = self._walberla_stencil[direction][0]
        cy = self._walberla_stencil[direction][1]
        cz = self._walberla_stencil[direction][2]
        fluid_cell = "Cell(it.x(), it.y(), it.z())"
        boundary_cell = f"Cell(it.x() + {cx}, it.y() + {cy}, it.z() + {cz})"
        init_element = f"elementInitialiser({fluid_cell}, {boundary_cell}, blocks, *block)"
        init_list = [f"const {self._dtype} q = (({self._dtype}) {init_element});", "element.q = q;"]

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return f"std::function<{self._dtype}(const Cell &, const Cell &, {self._blocks} elementInitialiser; "


class OutflowAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object, target=Target.CPU, field_name='pdfs',
                 pdfs_data_type=None, zeroth_timestep=None):
        assert isinstance(boundary_object, ExtrapolationOutflow)
        self._stencil = boundary_object.stencil
        self._lb_method = boundary_object.lb_method
        self._normal_direction = boundary_object.normal_direction
        self._field_name = field_name
        self._target = target
        self._dtype = BasicType(boundary_object.data_type).c_name
        if pdfs_data_type is None:
            self._pdfs_data_type = "real_t"
        else:
            pdfs_data_type = BasicType(pdfs_data_type)
            self._pdfs_data_type = pdfs_data_type.c_name

        self._streaming_pattern = boundary_object.streaming_pattern
        if zeroth_timestep:
            self._zeroth_timestep = zeroth_timestep
        else:
            self._zeroth_timestep = Timestep.EVEN if is_inplace(self._streaming_pattern) else Timestep.BOTH
        super(OutflowAdditionalDataHandler, self).__init__(stencil=stencil)

        assert sum([a != 0 for a in self._normal_direction]) == 1, \
            "The outflow boundary is only implemented for straight walls at the moment."

    @property
    def constructor_argument_name(self):
        return f"{self._field_name}CPUID_" if self._target == Target.GPU else ""

    @property
    def constructor_arguments(self):
        return f", BlockDataID {self._field_name}CPUID_" if self._target == Target.GPU else ""

    @property
    def initialiser_list(self):
        return f"{self._field_name}CPUID({self._field_name}CPUID_)," if self._target == Target.GPU else ""

    @property
    def additional_field_data(self):
        identifier = "CPU" if self._target == Target.GPU else ""
        return f"auto {self._field_name} = block->getData< field::GhostLayerField<{self._pdfs_data_type}, " \
               f"{len(self._stencil)}> >({self._field_name}{identifier}ID); "

    def data_initialisation(self, direction_index):
        pdf_acc = AccessPdfValues(self._stencil,
                                  streaming_pattern=self._streaming_pattern,
                                  timestep=self._zeroth_timestep,
                                  streaming_dir='out')

        init_list = []
        for key, value in self.get_init_dict(pdf_acc, direction_index).items():
            init_list.append(f"element.{key} = {self._dtype}( {self._field_name}->get({value}) );")

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return f"BlockDataID {self._field_name}CPUID;"

    @property
    def stencil_info(self):
        stencil_info = []
        for i, d in enumerate(self._stencil):
            if any([a != 0 and b != 0 and a == b for a, b in zip(self._normal_direction, d)]):
                direction = d if self._dim == 3 else d + (0,)
                stencil_info.append((i, direction, ", ".join([str(e) for e in direction])))
        return stencil_info

    def get_init_dict(self, pdf_accessor, direction_index):
        """The Extrapolation Outflow boundary needs additional data. This function provides a list of all values
        which have to be initialised"""
        position = ["it.x()", "it.y()", "it.z()"]
        direction = self._stencil[direction_index]
        inv_dir = self._stencil.index(inverse_direction(direction))

        tangential_offset = tuple(offset - normal for offset, normal in zip(direction, self._normal_direction))

        result = {}
        pos = []
        offsets = numeric_offsets(pdf_accessor.accs[inv_dir])
        for p, o, t in zip(position, offsets, tangential_offset):
            pos.append(p + f" + cell_idx_c({str(o + t)})")
        if self._dim == 2:
            pos.append("0")
        pos.append(str(numeric_index(pdf_accessor.accs[inv_dir])[0]))
        result['pdf'] = ', '.join(pos)

        pos = []
        for p, o, t in zip(position, offsets, tangential_offset):
            pos.append(p + f" + cell_idx_c({str(o + t)})")
        if self._dim == 2:
            pos.append("0")
        pos.append(str(numeric_index(pdf_accessor.accs[inv_dir])[0]))
        result['pdf_nd'] = ', '.join(pos)

        return result


class DiffusionDirichletAdditionalDataHandler(AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, DiffusionDirichlet)
        self._boundary_object = boundary_object
        super(DiffusionDirichletAdditionalDataHandler, self).__init__(stencil=stencil)

    @property
    def constructor_arguments(self):
        return ", std::function<real_t(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& " \
               "diffusionCallback "

    @property
    def initialiser_list(self):
        return "elementInitaliser(diffusionCallback),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, *_):
        init_list = ["real_t InitialisatonAdditionalData = elementInitaliser(Cell(it.x(), it.y(), it.z()), "
                     "blocks, *block);", "element.concentration = InitialisatonAdditionalData;"]

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return "std::function<real_t(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> " \
               "elementInitaliser; "
