from pystencils.stencil import inverse_direction


class AdditionalDataHandler:
    """Base class that defines how to handle boundary conditions holding additional data."""

    def __init__(self, stencil):
        self._dim = len(stencil[0])

        # waLBerla is a 3D framework. Therefore, a zero for the z index has to be added if we work in 2D
        if self._dim == 2:
            self._walberla_stencil = ()
            for d in stencil:
                d = d + (0,)
                self._walberla_stencil = self._walberla_stencil + (d,)
        else:
            self._walberla_stencil = stencil

    @property
    def constructor_arguments(self):
        return ""

    @property
    def initialiser_list(self):
        return ""

    @property
    def additional_arguments_for_fill_function(self):
        return ""

    @property
    def additional_parameters_for_fill_function(self):
        return ""

    @property
    def additional_field_data(self):
        return ""

    def data_initialisation(self, direction_index):
        return ""

    @property
    def additional_member_variable(self):
        return ""

    @property
    def stencil_info(self):
        return [(i, d, ", ".join([str(e) for e in d])) for i, d in enumerate(self._walberla_stencil)]

    @property
    def inverse_directions(self):
        inv_dirs = []
        for direction in self._walberla_stencil:
            inv_dirs.append(self._walberla_stencil.index(inverse_direction(direction)))
        return inv_dirs
