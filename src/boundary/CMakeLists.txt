###################################################################################################
#
# Module boundary
#
###################################################################################################


add_library( boundary )
target_link_libraries( boundary PUBLIC core communication domain_decomposition field stencil )
target_sources( boundary
      PRIVATE
      BoundaryHandling.h
      Boundary.h
      BoundaryHandlingCollection.h
      Boundary.cpp
      BoundaryUID.h
      )
