
###################################################################################################
#
# Module physics engine (pe)
#
###################################################################################################

waLBerla_add_module( DEPENDS core blockforest domain_decomposition geometry stencil vtk  )
target_link_libraries(pe PRIVATE lodepng)
