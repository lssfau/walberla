###################################################################################################
#
# Module walberla::pde
#
###################################################################################################

add_library( walberla_pde )
add_library( walberla::pde ALIAS walberla_pde )
target_link_libraries( walberla_pde PUBLIC walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::stencil )
target_sources( walberla_pde
      PRIVATE
      ConditionalResidualNorm.h
      ResidualNormStencilField.h
      ResidualNorm.h
      )

add_subdirectory( sweeps )
add_subdirectory( boundary )
add_subdirectory( iterations )
                         
###################################################################################################
