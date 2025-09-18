#############################################################################################################################
#
# Tests for generated walberla::lbm module
#
#############################################################################################################################
waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )

if( WALBERLA_BUILD_WITH_CODEGEN )

waLBerla_generate_target_from_python(NAME ExampleGenerated
        FILE Example.py
        CODEGEN_CFG example_codegen
        OUT_FILES LBMStorageSpecification.h LBMStorageSpecification.cpp
        LBMSweepCollection.h LBMSweepCollection.cpp
        NoSlip.h NoSlip.cpp
        UBB.h UBB.cpp
        LBMBoundaryCollection.h
        Example_InfoHeader.h)
waLBerla_add_test_executable ( GeneratedLBMExample Example.cpp )
target_link_libraries ( GeneratedLBMExample ExampleGenerated walberla::blockforest walberla::field
        walberla::lbm_generated walberla::timeloop )

waLBerla_generate_target_from_python(NAME InterpolationNoSlipGenerated
        FILE InterpolationNoSlip.py
        CODEGEN_CFG interpolation_no_slip_codegen
        OUT_FILES InterpolationNoSlipStorageSpecification.h InterpolationNoSlipStorageSpecification.cpp
        InterpolationNoSlipSweepCollection.h InterpolationNoSlipSweepCollection.cpp
        NoSlip.h NoSlip.cpp
        NoSlipBouzidi.h NoSlipBouzidi.cpp
        NoSlipQuadraticBB.h NoSlipQuadraticBB.cpp
        UBB.h UBB.cpp
        InterpolationNoSlipBoundaryCollection.h
        InterpolationNoSlipHeader.h)

waLBerla_add_test_executable ( InterpolationNoSlip InterpolationNoSlip.cpp )
target_link_libraries ( InterpolationNoSlip InterpolationNoSlipGenerated walberla::blockforest walberla::core
        walberla::field walberla::geometry walberla::lbm_generated walberla::timeloop )
# waLBerla_execute_test( NAME InterpolationNoSlip1 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm -Parameters.distanceWall=0.1 )
# waLBerla_execute_test( NAME InterpolationNoSlip2 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm -Parameters.distanceWall=0.5 )
waLBerla_execute_test( NAME InterpolationNoSlip3 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm )

waLBerla_generate_target_from_python(NAME FreeSlipRefinementGenerated
        FILE FreeSlipRefinement.py
        CODEGEN_CFG free_slip_refinement_codegen
        OUT_FILES FreeSlipRefinementStorageSpecification.h FreeSlipRefinementStorageSpecification.cpp
        FreeSlipRefinementSweepCollection.h FreeSlipRefinementSweepCollection.cpp
        FreeSlip.h FreeSlip.cpp
        UBB.h UBB.cpp
        Outflow.h Outflow.cpp
        FreeSlipRefinementBoundaryCollection.h
        FreeSlipRefinementInfoHeader.h)
waLBerla_add_test_executable ( FreeSlipRefinement FreeSlipRefinement.cpp )
target_link_libraries ( FreeSlipRefinement FreeSlipRefinementGenerated walberla::blockforest walberla::field
        walberla::lbm_generated walberla::timeloop )

if( WALBERLA_DOUBLE_ACCURACY )
   waLBerla_add_test_executable ( LDC LDC.cpp )
   target_link_libraries ( LDC walberla::blockforest walberla::field walberla::lbm_generated walberla::timeloop )
endif()

endif()
