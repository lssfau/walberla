add_subdirectory( parameter_files )

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))
        set(coupling_method 1)
        set(config "couplingMethod_${coupling_method}")
        message(STATUS "Config value is: ${config}" )
        waLBerla_generate_target_from_python(NAME Validations_3D_${config}
                FILE Validations_3D.py
                CODEGEN_CFG ${config}
                OUT_FILES LBMFluidSweep.h LBMFluidSweep.${CODEGEN_FILE_SUFFIX}
                LBMConcentrationSweep.h LBMConcentrationSweep.${CODEGEN_FILE_SUFFIX}
                LBMFluidSplitSweep.h LBMFluidSplitSweep.${CODEGEN_FILE_SUFFIX}
                LBMConcentrationSplitSweep.h LBMConcentrationSplitSweep.${CODEGEN_FILE_SUFFIX}
                PackInfoFluid.h PackInfoFluid.${CODEGEN_FILE_SUFFIX}
                PackInfoConcentration.h PackInfoConcentration.${CODEGEN_FILE_SUFFIX}
                InitializeFluidDomain.h InitializeFluidDomain.${CODEGEN_FILE_SUFFIX}
                InitializeConcentrationDomain.h InitializeConcentrationDomain.${CODEGEN_FILE_SUFFIX}
                BC_Fluid_NoSlip.h BC_Fluid_NoSlip.${CODEGEN_FILE_SUFFIX}
                BC_Fluid_UBB.h BC_Fluid_UBB.${CODEGEN_FILE_SUFFIX}
                BC_Fluid_Density.h BC_Fluid_Density.${CODEGEN_FILE_SUFFIX}
                BC_Fluid_FreeSlip.h BC_Fluid_FreeSlip.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_NoSlip.h BC_Concentration_NoSlip.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_UBB.h BC_Concentration_UBB.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_Density.h BC_Concentration_Density.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_FreeSlip.h BC_Concentration_FreeSlip.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_Neumann.h BC_Concentration_Neumann.${CODEGEN_FILE_SUFFIX}
                BC_Fluid_Neumann.h BC_Fluid_Neumann.${CODEGEN_FILE_SUFFIX}
                FluidMacroSetter.h FluidMacroSetter.cpp
                FluidMacroGetter.h FluidMacroGetter.cpp
                ConcentrationMacroSetter.h ConcentrationMacroSetter.cpp
                ConcentrationMacroGetter.h ConcentrationMacroGetter.cpp
                GeneralInfoHeader.h
        )

        waLBerla_add_executable(NAME 3d_gaussian_paper FILES 3D_Validations.cpp ../../utilities/InitializerFunctions.cpp
                DEPENDS blockforest core field geometry gpu lbm lbm_mesapd_coupling mesa_pd sqlite vtk Validations_3D_${config})
        target_compile_definitions(3d_gaussian_paper PRIVATE Weighting=2)
    endif ()
endif ()
