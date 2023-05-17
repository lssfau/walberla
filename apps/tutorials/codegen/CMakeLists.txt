# Code Generation Tutorials

if( WALBERLA_BUILD_WITH_CODEGEN )

    #   Tutorial 1: Heat Equation with pystencils
    walberla_generate_target_from_python( NAME CodegenHeatEquationKernel
        FILE HeatEquationKernel.py
        OUT_FILES HeatEquationKernel.cpp HeatEquationKernel.h )
    
    walberla_add_executable ( NAME 01_CodegenHeatEquation 
                              FILES 01_CodegenHeatEquation.cpp
                              DEPENDS blockforest core field stencil timeloop vtk pde CodegenHeatEquationKernel )

    #   Tutorial 2: lbmpy Lattice Model Generation
    waLBerla_link_files_to_builddir( *.prm )
    
    walberla_generate_target_from_python( NAME 02_LBMLatticeModelGenerationPython
        FILE 02_LBMLatticeModelGeneration.py
        OUT_FILES   SRTLatticeModel.cpp SRTLatticeModel.h 
                    SRTPackInfo.cpp SRTPackInfo.h )

    walberla_add_executable ( NAME 02_LBMLatticeModelGenerationApp 
                    FILES 02_LBMLatticeModelGeneration.cpp
                    DEPENDS blockforest core domain_decomposition field geometry timeloop lbm stencil vtk 02_LBMLatticeModelGenerationPython )

    #   Tutorial 3: Advanced lbmpy Code Generation

    walberla_generate_target_from_python( NAME 03_AdvancedLBMCodegenPython
        FILE 03_AdvancedLBMCodegen.py
        OUT_FILES   CumulantMRTSweep.${CODEGEN_FILE_SUFFIX} CumulantMRTSweep.h
                    CumulantMRTPackInfo.${CODEGEN_FILE_SUFFIX} CumulantMRTPackInfo.h
                    InitialPDFsSetter.${CODEGEN_FILE_SUFFIX} InitialPDFsSetter.h
                    CumulantMRTNoSlip.${CODEGEN_FILE_SUFFIX} CumulantMRTNoSlip.h)

    if(WALBERLA_BUILD_WITH_GPU_SUPPORT)
        walberla_add_executable ( NAME 03_AdvancedLBMCodegenApp
                        FILES 03_AdvancedLBMCodegen.cpp
                        DEPENDS blockforest gpu core domain_decomposition field geometry timeloop lbm stencil vtk 03_AdvancedLBMCodegenPython )
    else()
        walberla_add_executable ( NAME 03_AdvancedLBMCodegenApp
                FILES 03_AdvancedLBMCodegen.cpp
                DEPENDS blockforest core domain_decomposition field geometry timeloop lbm stencil vtk 03_AdvancedLBMCodegenPython )
    endif()

endif()