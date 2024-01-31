#############################################################################################################################
#
# Tests for generated lbm module
#
#############################################################################################################################
waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )

waLBerla_generate_target_from_python(NAME ExampleGenerated
        FILE Example.py
        OUT_FILES LBMStorageSpecification.h LBMStorageSpecification.cpp
        LBMSweepCollection.h LBMSweepCollection.cpp
        NoSlip.h NoSlip.cpp
        UBB.h UBB.cpp
        LBMBoundaryCollection.h
        Example_InfoHeader.h)
waLBerla_compile_test( FILES Example.cpp DEPENDS ExampleGenerated blockforest field lbm_generated timeloop )

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
