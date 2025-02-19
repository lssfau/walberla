waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )


foreach(streaming_pattern pull push aa esotwist esopull esopush)
    foreach(stencil d3q19 d3q27)
        foreach (collision_setup srt trt mrt mrt-overrelax central central-overrelax cumulant cumulant-overrelax cumulant-K17 entropic smagorinsky qr)
	    # KBC methods only for D2Q9 and D3Q27 defined
	    if (${collision_setup} STREQUAL "entropic" AND ${stencil} STREQUAL "d3q19")
		    continue()
        endif (${collision_setup} STREQUAL "entropic" AND ${stencil} STREQUAL "d3q19")
        if (${collision_setup} STREQUAL "cumulant-K17" AND ${stencil} STREQUAL "d3q19")
            continue()
        endif (${collision_setup} STREQUAL "cumulant-K17" AND ${stencil} STREQUAL "d3q19")
            set(config ${stencil}_${streaming_pattern}_${collision_setup})
            waLBerla_generate_target_from_python(NAME UniformGridCPUGenerated_${config}
                    FILE UniformGridCPU.py
                    CODEGEN_CFG ${config}
                    OUT_FILES UniformGridCPUStorageSpecification.h UniformGridCPUStorageSpecification.cpp
                    UniformGridCPUSweepCollection.h UniformGridCPUSweepCollection.cpp
                    NoSlip.cpp NoSlip.h
                    UBB.cpp UBB.h
                    UniformGridCPUBoundaryCollection.h
                    UniformGridCPU_StreamOnlyKernel.cpp UniformGridCPU_StreamOnlyKernel.h
                    UniformGridCPU_InfoHeader.h
                    )


            waLBerla_add_executable(NAME UniformGridCPU_${config}
                    FILES UniformGridCPU.cpp
                  DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::geometry walberla::python_coupling walberla::timeloop walberla::vtk UniformGridCPUGenerated_${config} )

            # all configs are excluded from all except for pull d3q27.
            if (${streaming_pattern} STREQUAL "pull" AND ${stencil} STREQUAL "d3q27")
                set_target_properties( UniformGridCPUGenerated_${config} PROPERTIES EXCLUDE_FROM_ALL FALSE)
                set_target_properties( UniformGridCPU_${config} PROPERTIES EXCLUDE_FROM_ALL FALSE)
            else()
                set_target_properties( UniformGridCPUGenerated_${config} PROPERTIES EXCLUDE_FROM_ALL TRUE)
                set_target_properties( UniformGridCPU_${config} PROPERTIES EXCLUDE_FROM_ALL TRUE)
            endif(${streaming_pattern} STREQUAL "pull" AND ${stencil} STREQUAL "d3q27")

        endforeach ()
    endforeach()
endforeach()
