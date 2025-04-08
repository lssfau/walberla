waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable ( NAME PE_GranularGas
                          FILES PE_GranularGas.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp
      DEPENDS walberla::blockforest walberla::core walberla::pe walberla::postprocessing walberla::sqlite )

waLBerla_add_executable ( NAME PE_LoadBalancing
                          FILES PE_LoadBalancing.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp
      DEPENDS walberla::blockforest walberla::core walberla::pe walberla::postprocessing walberla::sqlite )

waLBerla_add_executable ( NAME MESA_PD_LoadBalancing
                          FILES MESA_PD_LoadBalancing.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp sortParticleStorage.cpp CreateParticles.cpp
      DEPENDS walberla::blockforest walberla::core walberla::pe walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )

waLBerla_add_executable ( NAME MESA_PD_GranularGas
                          FILES MESA_PD_GranularGas.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp sortParticleStorage.cpp CreateParticles.cpp
      DEPENDS walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )

waLBerla_add_executable ( NAME MESA_PD_KernelBenchmark
                          FILES MESA_PD_KernelBenchmark.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp sortParticleStorage.cpp CreateParticles.cpp
      DEPENDS walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )

waLBerla_add_executable ( NAME MESA_PD_KernelLoadBalancing
                          FILES MESA_PD_KernelLoadBalancing.cpp SQLProperties.cpp Parameters.cpp NodeTimings.cpp sortParticleStorage.cpp CreateParticles.cpp
      DEPENDS walberla::blockforest walberla::core walberla::pe walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )
