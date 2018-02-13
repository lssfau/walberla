###################################################################################################
#
# Module mesh
#
###################################################################################################

waLBerla_add_module( DEPENDS blockforest boundary core domain_decomposition
                             python_coupling field geometry pe stencil vtk BUILD_ONLY_IF_FOUND OpenMesh )

###################################################################################################
