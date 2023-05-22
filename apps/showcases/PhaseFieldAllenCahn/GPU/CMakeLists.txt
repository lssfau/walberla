waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)
waLBerla_link_files_to_builddir(*.obj)

waLBerla_generate_target_from_python(NAME PhaseFieldCodeGenGPU
        FILE multiphase_codegen.py
        OUT_FILES initialize_phase_field_distributions.${CODEGEN_FILE_SUFFIX} initialize_phase_field_distributions.h
        initialize_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} initialize_velocity_based_distributions.h
        phase_field_LB_step.${CODEGEN_FILE_SUFFIX} phase_field_LB_step.h
        phase_field_LB_NoSlip.${CODEGEN_FILE_SUFFIX} phase_field_LB_NoSlip.h
        hydro_LB_step.${CODEGEN_FILE_SUFFIX} hydro_LB_step.h
        hydro_LB_NoSlip.${CODEGEN_FILE_SUFFIX} hydro_LB_NoSlip.h
        PackInfo_phase_field_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field_distributions.h
        PackInfo_phase_field.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field.h
        PackInfo_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_velocity_based_distributions.h
        ContactAngle.${CODEGEN_FILE_SUFFIX} ContactAngle.h
        GenDefines.h)

waLBerla_add_executable(NAME multiphaseGPU
        FILES multiphase.cpp PythonExports.cpp InitializerFunctions.cpp util.cpp multiphase_codegen.py
        DEPENDS blockforest core gpu field postprocessing python_coupling lbm geometry timeloop PhaseFieldCodeGenGPU)
