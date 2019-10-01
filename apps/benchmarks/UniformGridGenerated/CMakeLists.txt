waLBerla_python_file_generates(UniformGridGenerated.py
        UniformGridGenerated_LatticeModel.cpp
        UniformGridGenerated_Defines.h)


foreach(config trt )
    waLBerla_add_executable ( NAME UniformGridBenchmarkGenerated_${config}
            FILES UniformGridGenerated.cpp UniformGridGenerated.py
            DEPENDS blockforest boundary core domain_decomposition field geometry timeloop vtk gui
            CODEGEN_CFG ${config})
endforeach()

