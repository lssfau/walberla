
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

add_subdirectory( amr )
add_subdirectory( momentum_exchange_method )
add_subdirectory( utility )
add_subdirectory( mapping )