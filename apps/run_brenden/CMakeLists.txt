waLBerla_link_files_to_builddir( "run*" )

waLBerla_link_files_to_builddir( "*.conf" )

waLBerla_link_files_to_builddir( "*.py" )

#if ( WALBERLA_BUILD_WITH_OPENMESH )
   
    if( WALBERLA_BUILD_WITH_CODEGEN)

        waLBerla_generate_target_from_python(
                        NAME 	LBM_validation_py
                        FILE 	ScalarLBM.py
                                CODEGEN_CFG lbm_validation_codegen
                                OUT_FILES
                                PackInfo.h PackInfo.cpp
                                Neumann.h Neumann.cpp
                                SIConverter.h SIConverter.cpp
                                AdvectionDiffusionStorageSpecification.h AdvectionDiffusionStorageSpecification.cpp
                                AdvectionDiffusionSweepCollection.h AdvectionDiffusionSweepCollection.cpp
                                AdvectionDiffusionBoundaryCollection.h
                                InfoHeader.h

                        )	

        waLBerla_add_executable(
                        NAME 		LBM_validation
                        FILES		ScalarLBM.cpp
                        DEPENDS 	LBM_validation_py boundary core lbm lbm_generated  vtk
                )
        
    endif()
    #endif()
