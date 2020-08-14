# Code Generation Tutorials

if( WALBERLA_BUILD_WITH_CODEGEN )
    walberla_generate_target_from_python( NAME CodegenHeatEquationKernel
        FILE HeatEquationKernel.py
        OUT_FILES HeatEquationKernel.cpp HeatEquationKernel.h )
    
    walberla_add_executable ( NAME 01_CodegenHeatEquation 
                              FILES 01_CodegenHeatEquation.cpp
                              DEPENDS blockforest core field stencil timeloop vtk pde CodegenHeatEquationKernel )
endif()