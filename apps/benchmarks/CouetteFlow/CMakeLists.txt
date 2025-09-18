
waLBerla_link_files_to_builddir( "*.dat" )                 

add_executable ( CouetteFlow CouetteFlow.cpp)
target_link_libraries ( CouetteFlow walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm walberla::postprocessing walberla::stencil walberla::timeloop walberla::vtk walberla::sqlite )

##############
# Some tests #
##############

waLBerla_execute_test( NO_MODULE_LABEL NAME CouetteFlowTestNoCheckRelease COMMAND $<TARGET_FILE:CouetteFlow> ${CMAKE_CURRENT_SOURCE_DIR}/TestNoCheck.dat --trt --linear-exp PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo DEPENDS_ON_TARGETS CouetteFlow )
waLBerla_execute_test( NO_MODULE_LABEL NAME CouetteFlowTestNoCheckDebug   COMMAND $<TARGET_FILE:CouetteFlow> ${CMAKE_CURRENT_SOURCE_DIR}/TestNoCheck.dat --trt --linear-exp PROCESSES 4 LABELS longrun CONFIGURATIONS Debug DebugOptimized DEPENDS_ON_TARGETS CouetteFlow )

waLBerla_execute_test( NO_MODULE_LABEL NAME CouetteFlowTest0 COMMAND $<TARGET_FILE:CouetteFlow> ${CMAKE_CURRENT_SOURCE_DIR}/Test0.dat --trt --linear-exp LABELS longrun CONFIGURATIONS Release RelWithDbgInfo DEPENDS_ON_TARGETS CouetteFlow )
waLBerla_execute_test( NO_MODULE_LABEL NAME CouetteFlowTest2 COMMAND $<TARGET_FILE:CouetteFlow> ${CMAKE_CURRENT_SOURCE_DIR}/Test2.dat --trt --linear-exp LABELS longrun verylongrun PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo DEPENDS_ON_TARGETS CouetteFlow )
