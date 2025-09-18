add_executable ( PackPerformance PackPerformance.cpp )
target_link_libraries( PackPerformance walberla::core )

add_executable ( ProbeVsExtraMessage ProbeVsExtraMessage.cpp )
target_link_libraries ( ProbeVsExtraMessage walberla::core walberla::postprocessing walberla::stencil walberla::sqlite )
