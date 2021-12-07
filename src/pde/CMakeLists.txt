###################################################################################################
#
# Module pde 
#
###################################################################################################

add_library( pde )
target_link_libraries( pde PUBLIC blockforest boundary core domain_decomposition field stencil )
target_sources( pde
      PRIVATE
      ConditionalResidualNorm.h
      ResidualNormStencilField.h
      ResidualNorm.h
      )

add_subdirectory( sweeps )
add_subdirectory( boundary )
add_subdirectory( iterations )
                         
###################################################################################################
