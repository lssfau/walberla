###################################################################################################
#
# Module cuda
#
###################################################################################################

waLBerla_add_module( DEPENDS blockforest core communication domain_decomposition executiontree python_coupling field stencil
                     BUILD_ONLY_IF_FOUND CUDA )

###################################################################################################