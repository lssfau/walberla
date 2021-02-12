waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)
waLBerla_link_files_to_builddir(*.obj)

waLBerla_generate_target_from_python(NAME PhaseFieldCodeGenGPU
        FILE multiphase_codegen.py
        OUT_FILES initialize_phase_field_distributions.cu initialize_phase_field_distributions.h
        initialize_velocity_based_distributions.cu initialize_velocity_based_distributions.h
        phase_field_LB_step.cu phase_field_LB_step.h
        phase_field_LB_NoSlip.cu phase_field_LB_NoSlip.h
        hydro_LB_step.cu hydro_LB_step.h
        hydro_LB_NoSlip.cu hydro_LB_NoSlip.h
        stream_hydro.cu stream_hydro.h
        PackInfo_phase_field_distributions.cu PackInfo_phase_field_distributions.h
        PackInfo_phase_field.cu PackInfo_phase_field.h
        PackInfo_velocity_based_distributions.cu PackInfo_velocity_based_distributions.h
        GenDefines.h)

waLBerla_add_executable(NAME multiphaseGPU
        FILES multiphase.cpp PythonExports.cpp InitializerFunctions.cpp CalculateNormals.cpp contact.cu multiphase_codegen.py
        DEPENDS blockforest core cuda field postprocessing lbm geometry timeloop gui PhaseFieldCodeGenGPU)

set_target_properties(multiphaseGPU PROPERTIES CXX_VISIBILITY_PRESET hidden)
