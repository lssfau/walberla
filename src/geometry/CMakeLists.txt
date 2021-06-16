###################################################################################################
#
# Module geometry
#
###################################################################################################

waLBerla_add_module( DEPENDS core
                             domain_decomposition
                             boundary
                             field
                             python_coupling
                             stencil )
target_link_libraries(geometry PRIVATE lodepng)

###################################################################################################
