waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )


foreach(streaming_pattern pull push aa esotwist)
    foreach(stencil d3q19 d3q27)
        foreach (collision_setup srt trt mrt mrt-overrelax central central-overrelax cumulant cumulant-overrelax entropic smagorinsky)
	    # KBC methods only for D2Q9 and D3Q27 defined
	    if (${collision_setup} STREQUAL "entropic" AND ${stencil} STREQUAL "d3q19")
		    continue()
	    endif (${collision_setup} STREQUAL "entropic" AND ${stencil} STREQUAL "d3q19")
            set(config ${stencil}_${streaming_pattern}_${collision_setup})
            waLBerla_generate_target_from_python(NAME UniformGridGPUGenerated_${config}
                    FILE UniformGridGPU.py
                    CODEGEN_CFG ${config}
                    OUT_FILES   UniformGridGPU_LbKernel.${CODEGEN_FILE_SUFFIX} UniformGridGPU_LbKernel.h
                    UniformGridGPU_PackInfoEven.${CODEGEN_FILE_SUFFIX} UniformGridGPU_PackInfoEven.h
                    UniformGridGPU_PackInfoOdd.${CODEGEN_FILE_SUFFIX} UniformGridGPU_PackInfoOdd.h
                    UniformGridGPU_NoSlip.${CODEGEN_FILE_SUFFIX} UniformGridGPU_NoSlip.h
                    UniformGridGPU_UBB.${CODEGEN_FILE_SUFFIX} UniformGridGPU_UBB.h
                    UniformGridGPU_MacroSetter.${CODEGEN_FILE_SUFFIX} UniformGridGPU_MacroSetter.h
                    UniformGridGPU_StreamOnlyKernel.${CODEGEN_FILE_SUFFIX} UniformGridGPU_StreamOnlyKernel.h
                    UniformGridGPU_InfoHeader.h
                    )

            waLBerla_add_executable(NAME UniformGridGPU_${config}
                    FILES UniformGridGPU.cpp
                    DEPENDS blockforest boundary core gpu domain_decomposition field geometry python_coupling timeloop vtk UniformGridGPUGenerated_${config})

            # all configs are excluded from all except for pull d3q27.
            if (${streaming_pattern} STREQUAL "pull" AND ${stencil} STREQUAL "d3q27")
                set_target_properties( UniformGridGPUGenerated_${config} PROPERTIES EXCLUDE_FROM_ALL FALSE)
                set_target_properties( UniformGridGPU_${config} PROPERTIES EXCLUDE_FROM_ALL FALSE)
            else()
                set_target_properties( UniformGridGPUGenerated_${config} PROPERTIES EXCLUDE_FROM_ALL TRUE)
                set_target_properties( UniformGridGPU_${config} PROPERTIES EXCLUDE_FROM_ALL TRUE)
            endif(${streaming_pattern} STREQUAL "pull" AND ${stencil} STREQUAL "d3q27")

        endforeach ()
    endforeach()
endforeach()
