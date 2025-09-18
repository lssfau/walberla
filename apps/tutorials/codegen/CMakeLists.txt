# Code Generation Tutorials

if( WALBERLA_BUILD_WITH_CODEGEN )

    #   Tutorial 1: Heat Equation with pystencils
    walberla_generate_target_from_python( NAME CodegenHeatEquationKernel
        FILE HeatEquationKernel.py
        OUT_FILES HeatEquationKernel.cpp HeatEquationKernel.h )

    add_executable( 01_CodegenHeatEquation 01_CodegenHeatEquation.cpp )
    target_link_libraries ( 01_CodegenHeatEquation walberla::blockforest walberla::core walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde CodegenHeatEquationKernel )

    #   Tutorial 2: lbmpy Lattice Model Generation
    waLBerla_link_files_to_builddir( *.prm )
    
    walberla_generate_target_from_python( NAME 02_LBMLatticeModelGenerationPython
        FILE 02_LBMLatticeModelGeneration.py
        OUT_FILES   SRTLatticeModel.cpp SRTLatticeModel.h 
                    SRTPackInfo.cpp SRTPackInfo.h )

    add_executable ( 02_LBMLatticeModelGenerationApp 02_LBMLatticeModelGeneration.cpp )
    target_link_libraries( 02_LBMLatticeModelGenerationApp walberla::blockforest walberla::core
            walberla::domain_decomposition walberla::field walberla::geometry walberla::timeloop
            walberla::lbm walberla::stencil walberla::vtk 02_LBMLatticeModelGenerationPython )

    #   Tutorial 3: Advanced lbmpy Code Generation
    walberla_generate_target_from_python( NAME 03_AdvancedLBMCodegenPython
        FILE 03_AdvancedLBMCodegen.py
        OUT_FILES   CumulantMRTSweep.${CODEGEN_FILE_SUFFIX} CumulantMRTSweep.h
                    CumulantMRTPackInfo.${CODEGEN_FILE_SUFFIX} CumulantMRTPackInfo.h
                    InitialPDFsSetter.${CODEGEN_FILE_SUFFIX} InitialPDFsSetter.h
                    CumulantMRTNoSlip.${CODEGEN_FILE_SUFFIX} CumulantMRTNoSlip.h)

    if(WALBERLA_BUILD_WITH_GPU_SUPPORT)
        add_executable ( 03_AdvancedLBMCodegenApp 03_AdvancedLBMCodegen.cpp )
        target_link_libraries ( 03_AdvancedLBMCodegenApp walberla::blockforest walberla::gpu walberla::core
                walberla::domain_decomposition walberla::field walberla::geometry walberla::timeloop walberla::lbm
                walberla::stencil walberla::vtk 03_AdvancedLBMCodegenPython )
    else()
        add_executable ( 03_AdvancedLBMCodegenApp 03_AdvancedLBMCodegen.cpp )
        target_link_libraries ( 03_AdvancedLBMCodegenApp walberla::blockforest walberla::core
                walberla::domain_decomposition walberla::field walberla::geometry
                walberla::timeloop walberla::lbm walberla::stencil walberla::vtk 03_AdvancedLBMCodegenPython )
    endif()

    #   Tutorial 4: Advanced lbmpy Code Generation with Complex geometry
    if ( WALBERLA_BUILD_WITH_OPENMESH )

        # Links obj files in PROJECT_ROOT_DIR/geometry to build dir.
        waLBerla_link_geometry_to_builddir(cube.obj)
        waLBerla_link_geometry_to_builddir(bunny.obj)
        
        waLBerla_link_files_to_builddir( "04_LBComplexGeometryConfig.py" )

        set(TUT4_OUT_FILES)
        list(APPEND TUT4_OUT_FILES
            FreeSlipBC.h FreeSlipBC.${CODEGEN_FILE_SUFFIX}
            NoSlipBC.h NoSlipBC.${CODEGEN_FILE_SUFFIX}
            NoSlipQBBBC.h NoSlipQBBBC.${CODEGEN_FILE_SUFFIX}
            ObjNoSlipQBBBC.h ObjNoSlipQBBBC.${CODEGEN_FILE_SUFFIX}
            UBBBC.h UBBBC.${CODEGEN_FILE_SUFFIX}
            OutflowBC.h OutflowBC.${CODEGEN_FILE_SUFFIX}
            LBComplexGeometryStorageSpecification.h LBComplexGeometryStorageSpecification.${CODEGEN_FILE_SUFFIX}
            LBComplexGeometrySweepCollection.h LBComplexGeometrySweepCollection.${CODEGEN_FILE_SUFFIX}
            LBComplexGeometryBoundaryCollection.h
            InfoHeader.h
        )

        set( TUT4_DEPENDENCIES walberla::blockforest walberla::boundary walberla::core walberla::field
                walberla::geometry walberla::lbm walberla::lbm_generated walberla::mesh walberla::python_coupling
                walberla::stencil walberla::timeloop walberla::vtk 04_LBComplexGeometryPython )
        
        if(WALBERLA_BUILD_WITH_GPU_SUPPORT)
            list(APPEND TUT4_DEPENDENCIES walberla::gpu)
        endif()
        
        walberla_generate_target_from_python( 
                NAME 04_LBComplexGeometryPython
                FILE 04_LBComplexGeometry.py
                CODEGEN_CFG 04_codegen
                OUT_FILES  ${TUT4_OUT_FILES}
            )

        add_executable ( 04_LBComplexGeometryCodeGen 04_LBComplexGeometry.cpp )
        target_link_libraries( 04_LBComplexGeometryCodeGen ${TUT4_DEPENDENCIES} )
    endif()

endif()