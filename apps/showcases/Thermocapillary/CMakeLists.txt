waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)
waLBerla_link_files_to_builddir(*.obj)

waLBerla_generate_target_from_python(NAME ThermocapillaryGen
        FILE thermocapillary_codegen.py
        OUT_FILES initialize_phase_field_distributions.${CODEGEN_FILE_SUFFIX} initialize_phase_field_distributions.h
        initialize_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} initialize_velocity_based_distributions.h
        initialize_thermal_distributions.${CODEGEN_FILE_SUFFIX} initialize_thermal_distributions.h
        phase_field_LB_step.${CODEGEN_FILE_SUFFIX} phase_field_LB_step.h
        phase_field_LB_NoSlip.${CODEGEN_FILE_SUFFIX} phase_field_LB_NoSlip.h
        hydro_LB_step.${CODEGEN_FILE_SUFFIX} hydro_LB_step.h
        hydro_LB_NoSlip.${CODEGEN_FILE_SUFFIX} hydro_LB_NoSlip.h
        hydro_LB_UBB.${CODEGEN_FILE_SUFFIX} hydro_LB_UBB.h
        thermal_LB_step.${CODEGEN_FILE_SUFFIX} thermal_LB_step.h
        thermal_LB_DiffusionDirichlet_static.${CODEGEN_FILE_SUFFIX} thermal_LB_DiffusionDirichlet_static.h
        thermal_LB_DiffusionDirichlet_dynamic.${CODEGEN_FILE_SUFFIX} thermal_LB_DiffusionDirichlet_dynamic.h
        thermal_LB_NeumannByCopy.${CODEGEN_FILE_SUFFIX} thermal_LB_NeumannByCopy.h
        initialize_rk2.${CODEGEN_FILE_SUFFIX} initialize_rk2.h
        rk2_first_stage.${CODEGEN_FILE_SUFFIX} rk2_first_stage.h
        rk2_second_stage.${CODEGEN_FILE_SUFFIX} rk2_second_stage.h
        initialize_rk4.${CODEGEN_FILE_SUFFIX} initialize_rk4.h
        rk4_first_stage.${CODEGEN_FILE_SUFFIX} rk4_first_stage.h
        rk4_second_stage.${CODEGEN_FILE_SUFFIX} rk4_second_stage.h
        rk4_third_stage.${CODEGEN_FILE_SUFFIX} rk4_third_stage.h
        rk4_fourth_stage.${CODEGEN_FILE_SUFFIX} rk4_fourth_stage.h
        PackInfo_phase_field_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field_distributions.h
        PackInfo_velocity_based_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_velocity_based_distributions.h
        PackInfo_thermal_field_distributions.${CODEGEN_FILE_SUFFIX} PackInfo_thermal_field_distributions.h
        PackInfo_phase_field.${CODEGEN_FILE_SUFFIX} PackInfo_phase_field.h
        PackInfo_temperature_field.${CODEGEN_FILE_SUFFIX} PackInfo_temperature_field.h
        ContactAngle.${CODEGEN_FILE_SUFFIX} ContactAngle.h
        GenDefines.h)

if ( WALBERLA_BUILD_WITH_GPU_SUPPORT )
    waLBerla_add_executable(NAME thermocapillary
            FILES thermocapillary.cpp PythonExports.cpp InitializerFunctions.cpp thermocapillary_codegen.py
            DEPENDS blockforest core gpu field postprocessing python_coupling lbm geometry timeloop ThermocapillaryGen)
else ()
    waLBerla_add_executable(NAME thermocapillary
            FILES thermocapillary.cpp PythonExports.cpp InitializerFunctions.cpp thermocapillary_codegen.py
            DEPENDS blockforest core field postprocessing python_coupling lbm geometry timeloop ThermocapillaryGen)
endif ( WALBERLA_BUILD_WITH_GPU_SUPPORT )
