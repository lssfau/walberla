###################################################################################################
#
# Module lbm 
#
###################################################################################################

waLBerla_add_module( DEPENDS blockforest
                             boundary
                             communication
                             core
                             domain_decomposition
                             field
                             geometry
                             python_coupling
                             gui
                             stencil
                             timeloop
                             vtk )
                         
###################################################################################################

 