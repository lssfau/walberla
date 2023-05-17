 add_library( python_coupling )
 target_link_libraries( python_coupling PUBLIC pybind11::embed core communication domain_decomposition stencil field blockforest vtk )
 if( WALBERLA_BUILD_WITH_GPU_SUPPORT )
  target_link_libraries( python_coupling PUBLIC gpu )
 endif()
 target_sources( python_coupling
       PRIVATE
       Manager.cpp
       PythonWrapper.h
       DictWrapper.impl.h
       DictWrapper.h
       CreateConfig.h
       PythonCallback.cpp
       Manager.h
       PythonCallback.h
       CreateConfig.cpp
       )

 add_subdirectory( export )
 add_subdirectory( helper )