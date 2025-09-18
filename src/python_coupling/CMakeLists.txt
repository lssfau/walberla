add_library( walberla_python_coupling )
add_library( walberla::python_coupling ALIAS walberla_python_coupling )
target_link_libraries( walberla_python_coupling PUBLIC pybind11::embed walberla::core walberla::communication
        walberla::domain_decomposition walberla::stencil walberla::field walberla::blockforest walberla::vtk )

if ( WALBERLA_BUILD_WITH_GPU_SUPPORT )
    target_link_libraries( walberla_python_coupling PUBLIC walberla::gpu )
endif ()

target_sources( walberla_python_coupling
        PUBLIC
        PythonWrapper.h
        DictWrapper.impl.h
        DictWrapper.h
        CreateConfig.h
        Manager.h
        PythonCallback.h
)

target_sources( walberla_python_coupling
        PRIVATE
        Manager.cpp
        PythonCallback.cpp
        CreateConfig.cpp
)

add_subdirectory( export )
add_subdirectory( helper )