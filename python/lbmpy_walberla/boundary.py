import pystencils_walberla.boundary
from lbmpy.boundaries.boundaryhandling import create_lattice_boltzmann_boundary_kernel


def generate_boundary(generation_context,
                      class_name,
                      boundary_object,
                      lb_method,
                      field_name='pdfs',
                      **create_kernel_params):

    def boundary_creation_function(field, index_field, stencil, boundary_functor, target='cpu', openmp=True, **kwargs):
        return create_lattice_boltzmann_boundary_kernel(field,
                                                        index_field,
                                                        lb_method,
                                                        boundary_functor,
                                                        target=target,
                                                        **kwargs)

    pystencils_walberla.boundary.generate_boundary(generation_context,
                                                   class_name,
                                                   boundary_object,
                                                   field_name=field_name,
                                                   neighbor_stencil=lb_method.stencil,
                                                   index_shape=[len(lb_method.stencil)],
                                                   kernel_creation_function=boundary_creation_function,
                                                   namespace='lbm',
                                                   **create_kernel_params)
