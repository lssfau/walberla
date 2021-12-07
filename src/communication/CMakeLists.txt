###################################################################################################
#
# Module communication
#
###################################################################################################

add_library( communication INTERFACE )
target_sources( communication
      INTERFACE
      ReducePackInfo.h
      UniformMPIDatatypeInfo.h
      UniformPackInfo.h
      )
target_link_libraries( communication INTERFACE core domain_decomposition stencil )

###################################################################################################                        