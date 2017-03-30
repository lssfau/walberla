
###################################################################################################
#
# Module field ( lattice data structure )
#
###################################################################################################

OPTION ( WALBERLA_FIELD_MONITORED_ACCESS "Enable (slow) monitoring mechanism for Fields" )
OPTION ( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION "Enables/Disables thread-safe field allocation (do NOT disable if you are not absolutely sure what you are doing!)" ON )
mark_as_advanced( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION )

configure_file ( CMakeDefs.in.h  CMakeDefs.h )

waLBerla_add_module( DEPENDS blockforest core communication domain_decomposition python_coupling stencil vtk )

 