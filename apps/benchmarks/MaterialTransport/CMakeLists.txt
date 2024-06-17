waLBerla_link_files_to_builddir("*.prm")

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))
        waLBerla_add_executable(NAME MaterialTransport FILES MaterialTransport.cpp
                DEPENDS blockforest core field geometry gpu lbm lbm_mesapd_coupling mesa_pd sqlite vtk PSMCodegenPython_trt-smagorinsky_sc1)
        target_compile_definitions(MaterialTransport PRIVATE Weighting=2)
    endif ()
endif ()
