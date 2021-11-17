if (WALBERLA_BUILD_WITH_PYTHON)
   waLBerla_add_module(DEPENDS core communication domain_decomposition stencil field blockforest vtk cuda)
   target_link_libraries(python_coupling PUBLIC pybind11::embed core domain_decomposition)
endif()