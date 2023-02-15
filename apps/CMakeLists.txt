
# Benchmarks
if ( WALBERLA_BUILD_BENCHMARKS )
    add_subdirectory ( benchmarks )
endif()


# Tools
if ( WALBERLA_BUILD_TOOLS )
    add_subdirectory ( tools )
endif()


# Tutorials
if ( WALBERLA_BUILD_TUTORIALS )
    add_subdirectory ( tutorials )
endif()

# Showcases
if ( WALBERLA_BUILD_SHOWCASES )
    add_subdirectory ( showcases )
endif()

# Python module
if ( WALBERLA_BUILD_WITH_PYTHON )
    add_subdirectory( pythonmodule )
endif()
