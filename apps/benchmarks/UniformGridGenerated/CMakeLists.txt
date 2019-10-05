waLBerla_link_files_to_builddir( "*.prm" )


waLBerla_python_file_generates(UniformGridGenerated.py
        GenMacroGetter.cpp GenMacroSetter.cpp
        GenPackInfo.cpp GenPackInfoAAPush.cpp GenPackInfoAAPull.cpp
        GenLbKernel.cpp GenLbKernelAAEven.cpp GenLbKernelAAOdd.cpp
        GenDefines.h)


foreach(config trt )
    waLBerla_add_executable ( NAME UniformGridBenchmarkGenerated_${config}
            FILES UniformGridGenerated.cpp UniformGridGenerated.py
            DEPENDS blockforest boundary core domain_decomposition field geometry timeloop vtk gui
            CODEGEN_CFG ${config})
endforeach()

