waLBerla_add_executable ( NAME PackPerformance 
                          FILES PackPerformance.cpp
                          DEPENDS core )

waLBerla_add_executable ( NAME ProbeVsExtraMessage 
                          FILES ProbeVsExtraMessage.cpp
                          DEPENDS core postprocessing stencil sqlite)
