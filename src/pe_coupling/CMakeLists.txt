###################################################################################################
#
# Module pe_coupling 
#
###################################################################################################

add_library( pe_coupling )
target_link_libraries( pe_coupling PUBLIC blockforest boundary communication core domain_decomposition field geometry lbm pe stencil vtk )

add_subdirectory( amr )
add_subdirectory( discrete_particle_methods )
add_subdirectory( momentum_exchange_method )
add_subdirectory( utility )
add_subdirectory( partially_saturated_cells_method )
add_subdirectory( geometry )
add_subdirectory( mapping )
                         
###################################################################################################
