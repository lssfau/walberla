waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)

if (WALBERLA_BUILD_WITH_CUDA)
    waLBerla_generate_target_from_python(NAME BenchmarkPhaseFieldCodeGenGPU
            FILE multiphase_codegen.py
            OUT_FILES initialize_phase_field_distributions.cu initialize_phase_field_distributions.h
            initialize_velocity_based_distributions.cu initialize_velocity_based_distributions.h
            phase_field_LB_step.cu phase_field_LB_step.h
            hydro_LB_step.cu hydro_LB_step.h
            PackInfo_phase_field_distributions.cu PackInfo_phase_field_distributions.h
            PackInfo_phase_field.cu PackInfo_phase_field.h
            PackInfo_velocity_based_distributions.cu PackInfo_velocity_based_distributions.h
            GenDefines.h)

    waLBerla_add_executable(NAME benchmark_multiphase
            FILES benchmark_multiphase.cpp InitializerFunctions.cpp multiphase_codegen.py
            DEPENDS blockforest core cuda field postprocessing lbm geometry timeloop gui BenchmarkPhaseFieldCodeGenGPU)
    set_target_properties(benchmark_multiphase PROPERTIES CXX_VISIBILITY_PRESET hidden)
else ()
    waLBerla_generate_target_from_python(NAME BenchmarkPhaseFieldCodeGenCPU
            FILE multiphase_codegen.py
            OUT_FILES initialize_phase_field_distributions.cpp initialize_phase_field_distributions.h
            initialize_velocity_based_distributions.cpp initialize_velocity_based_distributions.h
            phase_field_LB_step.cpp phase_field_LB_step.h
            hydro_LB_step.cpp hydro_LB_step.h
            PackInfo_phase_field_distributions.cpp PackInfo_phase_field_distributions.h
            PackInfo_phase_field.cpp PackInfo_phase_field.h
            PackInfo_velocity_based_distributions.cpp PackInfo_velocity_based_distributions.h
            GenDefines.h)

    waLBerla_add_executable(NAME benchmark_multiphase
            FILES benchmark_multiphase.cpp InitializerFunctions.cpp multiphase_codegen.py
            DEPENDS blockforest core field postprocessing lbm geometry timeloop gui BenchmarkPhaseFieldCodeGenCPU)
    set_target_properties(benchmark_multiphase PROPERTIES CXX_VISIBILITY_PRESET hidden)
endif (WALBERLA_BUILD_WITH_CUDA)

