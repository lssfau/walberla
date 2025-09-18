waLBerla_link_files_to_builddir("*.prm")

if ( WALBERLA_BUILD_WITH_CODEGEN )
    if (
            NOT WALBERLA_BUILD_WITH_GPU_SUPPORT
            OR (
                WALBERLA_BUILD_WITH_GPU_SUPPORT
                AND ( CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP )
            )
    )
        add_executable ( Percolation Percolation.cpp )
        target_link_libraries ( Percolation walberla::blockforest walberla::core walberla::field
                walberla::geometry walberla::gpu walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
                walberla::sqlite walberla::vtk PSMCodegenPython_trt-smagorinsky_sc1 )
        target_compile_definitions ( Percolation PRIVATE Weighting=2 )
    endif ()
endif ()
