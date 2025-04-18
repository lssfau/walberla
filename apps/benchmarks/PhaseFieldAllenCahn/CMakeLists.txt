waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)

waLBerla_generate_target_from_python(NAME BenchmarkPhaseFieldCodeGen
        FILE multiphase_codegen.py
        OUT_FILES initialize_phase_field_distributions.${CODEGEN_FILE_SUFFIX} initialize_phase_field_distributions.h
        initialize_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} initialize_velocity_based_distributions.h
        phase_field_LB_step.${CODEGEN_FILE_SUFFIX} phase_field_LB_step.h
        hydro_LB_step.${CODEGEN_FILE_SUFFIX} hydro_LB_step.h
        PackInfo_phase_field_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field_distributions.h
        PackInfo_phase_field.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field.h
        PackInfo_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_velocity_based_distributions.h
        GenDefines.h)

if (WALBERLA_BUILD_WITH_GPU_SUPPORT )
    waLBerla_add_executable(NAME benchmark_multiphase
            FILES benchmark_multiphase.cpp InitializerFunctions.cpp multiphase_codegen.py
          DEPENDS walberla::blockforest walberla::core walberla::gpu walberla::field walberla::postprocessing walberla::python_coupling walberla::lbm_generated walberla::geometry walberla::timeloop BenchmarkPhaseFieldCodeGen )
else ()
    waLBerla_add_executable(NAME benchmark_multiphase
            FILES benchmark_multiphase.cpp InitializerFunctions.cpp multiphase_codegen.py
          DEPENDS walberla::blockforest walberla::core walberla::field walberla::postprocessing walberla::python_coupling walberla::lbm_generated walberla::geometry walberla::timeloop BenchmarkPhaseFieldCodeGen )
endif (WALBERLA_BUILD_WITH_GPU_SUPPORT )

