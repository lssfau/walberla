#############################################################################################################################
#
# Tests for generated lbm module
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
waLBerla_compile_test( FILES Example.cpp DEPENDS ExampleGenerated blockforest field lbm_generated timeloop )

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

waLBerla_compile_test( FILES InterpolationNoSlip.cpp DEPENDS InterpolationNoSlipGenerated blockforest core field geometry lbm_generated timeloop )
# waLBerla_execute_test( NAME InterpolationNoSlip1 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm -Parameters.distanceWall=0.1 )
# waLBerla_execute_test( NAME InterpolationNoSlip2 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm -Parameters.distanceWall=0.5 )
waLBerla_execute_test( NAME InterpolationNoSlip3 COMMAND $<TARGET_FILE:InterpolationNoSlip> ${CMAKE_CURRENT_SOURCE_DIR}/InterpolationNoSlip.prm )
endif()

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
waLBerla_compile_test( FILES FreeSlipRefinement.cpp DEPENDS FreeSlipRefinementGenerated blockforest field lbm_generated timeloop )

if( WALBERLA_DOUBLE_ACCURACY )
waLBerla_compile_test( FILES LDC.cpp DEPENDS blockforest field lbm_generated timeloop )
endif()
