import pystencils_walberla.boundary
from lbmpy.boundaries.boundaryhandling import create_lattice_boltzmann_boundary_kernel
from lbmpy.advanced_streaming import Timestep, is_inplace


def generate_boundary(generation_context,
                      class_name,
                      boundary_object,
                      lb_method,
                      field_name='pdfs',
                      streaming_pattern='pull',
                      always_generate_separate_classes=False,
                      additional_data_handler=None,
                      **create_kernel_params):
    def boundary_creation_function(field, index_field, stencil, boundary_functor, target='cpu', openmp=True, **kwargs):
        pargs = (field, index_field, lb_method, boundary_functor)
        kwargs = {'target': target, **kwargs}
        if is_inplace(streaming_pattern) or always_generate_separate_classes:
            return {
                'EvenSweep': create_lattice_boltzmann_boundary_kernel(*pargs,
                                                                      streaming_pattern=streaming_pattern,
                                                                      prev_timestep=Timestep.EVEN,
                                                                      **kwargs),
                'OddSweep': create_lattice_boltzmann_boundary_kernel(*pargs,
                                                                     streaming_pattern=streaming_pattern,
                                                                     prev_timestep=Timestep.ODD,
                                                                     **kwargs)
            }
        else:
            return create_lattice_boltzmann_boundary_kernel(*pargs,
                                                            streaming_pattern=streaming_pattern,
                                                            prev_timestep=Timestep.BOTH,
                                                            **kwargs)

    pystencils_walberla.boundary.generate_boundary(generation_context,
                                                   class_name,
                                                   boundary_object,
                                                   field_name=field_name,
                                                   neighbor_stencil=lb_method.stencil,
                                                   index_shape=[len(lb_method.stencil)],
                                                   kernel_creation_function=boundary_creation_function,
                                                   namespace='lbm',
                                                   additional_data_handler=additional_data_handler,
                                                   **create_kernel_params)
