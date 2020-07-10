if ( NOT WALBERLA_BUILD_WITH_CUDA)
add_subdirectory( CPU )
endif()

if ( WALBERLA_BUILD_WITH_CUDA)
add_subdirectory( GPU )
endif()
