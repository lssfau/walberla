waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )


foreach(config trt smagorinsky mrt entropic_kbc_n4 cumulant )
    waLBerla_generate_target_from_python(NAME UniformGridGenerated_${config}
          CODEGEN_CFG ${config}
          FILE UniformGridGenerated.py
          OUT_FILES GenMacroGetter.cpp GenMacroGetter.h GenMacroSetter.cpp GenMacroSetter.h
          GenPackInfo.cpp GenPackInfo.h GenPackInfoAAPush.cpp GenPackInfoAAPush.h GenPackInfoAAPull.cpp GenPackInfoAAPull.h
          GenLbKernel.cpp GenLbKernel.h GenLbKernelAAEven.cpp GenLbKernelAAEven.h GenLbKernelAAOdd.cpp GenLbKernelAAOdd.h
          GenMpiDtypeInfo.h GenMpiDtypeInfoAAPull.h GenMpiDtypeInfoAAPush.h
          GenDefines.h)

    waLBerla_add_executable ( NAME UniformGridBenchmarkGenerated_${config}
          FILES UniformGridGenerated.cpp
          DEPENDS blockforest boundary core domain_decomposition field geometry timeloop vtk gui
          UniformGridGenerated_${config})

endforeach()
