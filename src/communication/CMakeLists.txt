###################################################################################################
#
# Module communication
#
###################################################################################################

add_library( walberla_communication INTERFACE )
add_library( walberla::communication ALIAS walberla_communication )
target_sources( walberla_communication
      INTERFACE
      ReducePackInfo.h
      UniformMPIDatatypeInfo.h
      UniformPackInfo.h
      )
target_link_libraries( walberla_communication INTERFACE walberla::core walberla::domain_decomposition walberla::stencil )

###################################################################################################                        