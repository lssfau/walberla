add_subdirectory( CPU )

if ( WALBERLA_BUILD_WITH_CUDA)
add_subdirectory( GPU )
endif()
