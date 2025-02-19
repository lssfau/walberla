waLBerla_link_files_to_builddir("*.prm")

if (WALBERLA_BUILD_WITH_GPU_SUPPORT AND WALBERLA_BUILD_WITH_CODEGEN AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP))
    waLBerla_add_executable(NAME FluidizedBed_PSM_GPU FILES FluidizedBedGPU.cpp
          DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::gpu walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk PSMCodegenPython_srt_sc1 )
endif ()
