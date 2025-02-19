
###################################################################################################
#
# Module walberla::lbm - mesapd - coupling
#
###################################################################################################

add_library( walberla_lbm_mesapd_coupling )
add_library( walberla::lbm_mesapd_coupling ALIAS walberla_lbm_mesapd_coupling )
target_link_libraries( walberla_lbm_mesapd_coupling PUBLIC walberla::boundary walberla::blockforest walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::stencil )
target_sources( walberla_lbm_mesapd_coupling
      PRIVATE
      DataTypes.h
      )

# Maximum number of particles that may overlap with a cell. For fully resolved particles, 2 should normally be
# sufficient (for a sufficiently high stiffness in the DEM).
set(MAX_PARTICLES_PER_CELL 2)
add_custom_target(MAX_PARTICLES_PER_CELL) # Make it a target such that the code generation runs again if changed
target_sources( walberla_lbm_mesapd_coupling
        PRIVATE
        DataTypesCodegen.h
        )
target_compile_definitions( walberla_lbm_mesapd_coupling PUBLIC MAX_PARTICLES_PER_CELL=${MAX_PARTICLES_PER_CELL} )


add_subdirectory( amr )
add_subdirectory( momentum_exchange_method )
add_subdirectory( overlapping )
add_subdirectory( partially_saturated_cells_method )
add_subdirectory( utility )
add_subdirectory( mapping )
