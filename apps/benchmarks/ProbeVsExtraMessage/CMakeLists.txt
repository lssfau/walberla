waLBerla_add_executable ( NAME PackPerformance 
                          FILES PackPerformance.cpp
      DEPENDS walberla::core )

waLBerla_add_executable ( NAME ProbeVsExtraMessage 
                          FILES ProbeVsExtraMessage.cpp
      DEPENDS walberla::core walberla::postprocessing walberla::stencil walberla::sqlite )
