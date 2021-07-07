
waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )


foreach(streaming_pattern aa) # choose from {pull, push, aa, esotwist}
    foreach(stencil d3q27) # choose from {d3q19 d3q27}
        foreach (collision_setup srt trt mrt cumulant) # choose from {srt trt mrt cumulant entropic smagorinsky}
            set(config ${stencil}_${streaming_pattern}_${collision_setup})
            waLBerla_generate_target_from_python(NAME UniformGridGPUGenerated_${config}
                    FILE UniformGridGPU.py
                    CODEGEN_CFG ${config}
                    OUT_FILES   UniformGridGPU_LbKernel.cu UniformGridGPU_LbKernel.h
                    UniformGridGPU_PackInfoEven.cu UniformGridGPU_PackInfoEven.h
                    UniformGridGPU_PackInfoOdd.cu UniformGridGPU_PackInfoOdd.h
                    UniformGridGPU_NoSlip.cu UniformGridGPU_NoSlip.h
                    UniformGridGPU_UBB.cu UniformGridGPU_UBB.h
                    UniformGridGPU_MacroSetter.cu UniformGridGPU_MacroSetter.h
                    UniformGridGPU_InfoHeader.h
                    )


            waLBerla_add_executable(NAME UniformGridGPU_${config}
                    FILES UniformGridGPU.cpp
                    DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk UniformGridGPUGenerated_${config})
            set_target_properties( UniformGridGPU_${config} PROPERTIES CXX_VISIBILITY_PRESET hidden)
        endforeach ()
    endforeach()
endforeach()