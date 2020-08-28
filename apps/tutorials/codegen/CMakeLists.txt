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
    walberla_generate_target_from_python( NAME LBMLatticeModelGenerationPython
        FILE LBMLatticeModelGeneration.py
        OUT_FILES   SRTLatticeModel.cpp SRTLatticeModel.h 
                    SRTPackInfo.cpp SRTPackInfo.h
                    SRTNoSlip.cpp SRTNoSlip.h
                    CumulantMRTLatticeModel.cpp CumulantMRTLatticeModel.h 
                    CumulantMRTPackInfo.cpp CumulantMRTPackInfo.h
                    CumulantMRTNoSlip.cpp CumulantMRTNoSlip.h
                    EntropicMRTLatticeModel.cpp EntropicMRTLatticeModel.h
                    EntropicMRTPackInfo.cpp EntropicMRTPackInfo.h
                    EntropicMRTNoSlip.cpp EntropicMRTNoSlip.h )

    waLBerla_link_files_to_builddir( *.prm )

    walberla_add_executable ( NAME 02_LBMLatticeModelGenerationApp 
                    FILES 02_LBMLatticeModelGeneration.cpp
                    DEPENDS blockforest core domain_decomposition field geometry timeloop lbm stencil vtk LBMLatticeModelGenerationPython )

endif()