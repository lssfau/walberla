###################################################################################################
#
# Module walberla::pe_coupling
#
###################################################################################################

add_library( walberla_pe_coupling )
add_library( walberla::pe_coupling ALIAS walberla_pe_coupling )
target_link_libraries( walberla_pe_coupling PUBLIC walberla::blockforest walberla::boundary walberla::communication walberla::core walberla::domain_decomposition walberla::field walberla::geometry walberla::lbm walberla::pe walberla::stencil walberla::vtk )

add_subdirectory( amr )
add_subdirectory( discrete_particle_methods )
add_subdirectory( momentum_exchange_method )
add_subdirectory( utility )
add_subdirectory( partially_saturated_cells_method )
add_subdirectory( geometry )
add_subdirectory( mapping )
                         
###################################################################################################
