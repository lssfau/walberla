
###################################################################################################
#
# Module lbm - mesapd - coupling
#
###################################################################################################

add_library( lbm_mesapd_coupling )
target_link_libraries( lbm_mesapd_coupling PUBLIC boundary blockforest core domain_decomposition field lbm mesa_pd stencil )
target_sources( lbm_mesapd_coupling
      PRIVATE
      DataTypes.h
      )

# Maximum number of particles that may overlap with a cell. For fully resolved particles, 2 should normally be
# sufficient (for a sufficiently high stiffness in the DEM).
set(MAX_PARTICLES_PER_CELL 2)
add_custom_target(MAX_PARTICLES_PER_CELL) # Make it a target such that the code generation runs again if changed
target_sources( lbm_mesapd_coupling
        PRIVATE
        DataTypesCodegen.h
        )
target_compile_definitions(lbm_mesapd_coupling PUBLIC MAX_PARTICLES_PER_CELL=${MAX_PARTICLES_PER_CELL})


add_subdirectory( amr )
add_subdirectory( momentum_exchange_method )
add_subdirectory( overlapping )
add_subdirectory( partially_saturated_cells_method )
add_subdirectory( utility )
add_subdirectory( mapping )
