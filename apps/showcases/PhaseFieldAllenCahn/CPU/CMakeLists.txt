waLBerla_link_files_to_builddir(*.prm)
waLBerla_link_files_to_builddir(*.py)
waLBerla_link_files_to_builddir(*.obj)

waLBerla_generate_target_from_python(NAME PhaseFieldCodeGenCPU
        FILE multiphase_codegen.py
        OUT_FILES initialize_phase_field_distributions.cpp initialize_phase_field_distributions.h
        initialize_velocity_based_distributions.cpp initialize_velocity_based_distributions.h
        phase_field_LB_step.cpp phase_field_LB_step.h
        phase_field_LB_NoSlip.cpp phase_field_LB_NoSlip.h
        hydro_LB_step.cpp hydro_LB_step.h
        hydro_LB_NoSlip.cpp hydro_LB_NoSlip.h
        PackInfo_phase_field_distributions.cpp PackInfo_phase_field_distributions.h
        PackInfo_velocity_based_distributions.cpp PackInfo_velocity_based_distributions.h
        PackInfo_phase_field.cpp PackInfo_phase_field.h
        ContactAngle.cpp ContactAngle.h
        GenDefines.h)

waLBerla_add_executable(NAME multiphaseCPU
        FILES multiphase.cpp PythonExports.cpp InitializerFunctions.cpp multiphase_codegen.py
      DEPENDS walberla::blockforest walberla::core walberla::field walberla::postprocessing walberla::python_coupling walberla::lbm walberla::geometry walberla::timeloop PhaseFieldCodeGenCPU )
