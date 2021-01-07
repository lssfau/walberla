
waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )


foreach (config srt trt mrt smagorinsky entropic smagorinsky_noopt entropic_kbc_n4
      entropic_kbc_n4_noopt mrt_noopt mrt_full mrt_full_noopt
      cumulant cumulant_d3q27
      srt_d3q27 mrt_d3q27 mrt_d3q27_noopt smagorinsky_d3q27 smagorinsky_d3q27_noopt mrt_full_d3q27 mrt_full_d3q27_noopt)

    waLBerla_generate_target_from_python(NAME UniformGridGPUGenerated_${config}
          FILE UniformGridGPU.py
          CODEGEN_CFG ${config}
          OUT_FILES UniformGridGPU_LatticeModel.cpp UniformGridGPU_LatticeModel.h
          UniformGridGPU_LbKernel.cu UniformGridGPU_LbKernel.h
          UniformGridGPU_NoSlip.cu UniformGridGPU_NoSlip.h
          UniformGridGPU_UBB.cu UniformGridGPU_UBB.h
          UniformGridGPU_PackInfo.cu UniformGridGPU_PackInfo.h
          UniformGridGPU_MacroSetter.cpp UniformGridGPU_MacroSetter.h
          UniformGridGPU_MacroGetter.cpp UniformGridGPU_MacroGetter.h
          UniformGridGPU_Defines.h
          )


    waLBerla_add_executable(NAME UniformGridBenchmarkGPU_${config}
          FILES UniformGridGPU.cpp
          DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk gui UniformGridGPUGenerated_${config})
    set_target_properties( UniformGridBenchmarkGPU_${config} PROPERTIES CXX_VISIBILITY_PRESET hidden)
endforeach ()


foreach (config srt trt mrt smagorinsky entropic)

    waLBerla_generate_target_from_python(NAME UniformGridGPUGenerated_AA_${config}
          FILE UniformGridGPU_AA.py
          CODEGEN_CFG ${config}
          OUT_FILES UniformGridGPU_AA_PackInfoPull.cu UniformGridGPU_AA_PackInfoPull.h
          UniformGridGPU_AA_LbKernelOdd.cu UniformGridGPU_AA_LbKernelOdd.h
          UniformGridGPU_AA_LbKernelEven.cu UniformGridGPU_AA_LbKernelEven.h
          UniformGridGPU_AA_PackInfoPush.cu UniformGridGPU_AA_PackInfoPush.h
          UniformGridGPU_AA_MacroSetter.cpp UniformGridGPU_AA_MacroSetter.h
          UniformGridGPU_AA_MacroGetter.cpp UniformGridGPU_AA_MacroGetter.h
          UniformGridGPU_AA_Defines.h
          )


    waLBerla_add_executable(NAME UniformGridBenchmarkGPU_AA_${config}
          FILES UniformGridGPU_AA.cpp
          DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk gui UniformGridGPUGenerated_AA_${config})
    set_target_properties( UniformGridBenchmarkGPU_AA_${config} PROPERTIES CXX_VISIBILITY_PRESET hidden)
endforeach ()
