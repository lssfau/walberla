waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.py )

add_executable ( MESA_PD_LoadBalancing MESA_PD_LoadBalancing.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp
        sortParticleStorage.cpp CreateParticles.cpp )
target_link_libraries ( MESA_PD_LoadBalancing walberla::blockforest walberla::core walberla::mesa_pd
        walberla::postprocessing walberla::sqlite walberla::vtk )

add_executable ( MESA_PD_GranularGas MESA_PD_GranularGas.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp
        sortParticleStorage.cpp CreateParticles.cpp )
target_link_libraries ( MESA_PD_GranularGas walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )

add_executable ( MESA_PD_KernelBenchmark MESA_PD_KernelBenchmark.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp
        sortParticleStorage.cpp CreateParticles.cpp )
target_link_libraries ( MESA_PD_KernelBenchmark walberla::blockforest walberla::core walberla::mesa_pd
        walberla::postprocessing walberla::sqlite walberla::vtk )

add_executable ( MESA_PD_KernelLoadBalancing MESA_PD_KernelLoadBalancing.cpp SQLProperties.cpp Parameters.cpp
        NodeTimings.cpp sortParticleStorage.cpp CreateParticles.cpp )
target_link_libraries ( MESA_PD_KernelLoadBalancing walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )
