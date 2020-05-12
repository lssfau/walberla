###################################################################################################
#
# Module mesh
#
###################################################################################################

waLBerla_add_module( DEPENDS blockforest boundary core domain_decomposition
                             python_coupling field geometry pe stencil mesa_pd mesh_common
                     BUILD_ONLY_IF_FOUND OpenMesh)

###################################################################################################
