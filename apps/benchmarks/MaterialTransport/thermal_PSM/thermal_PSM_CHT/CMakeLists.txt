waLBerla_link_files_to_builddir("*.prm")

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))
        foreach(solid_collision 1 2)

            set(config sc_${solid_collision})
            message(STATUS "Config value in thermal PSM is: ${config}" )
            waLBerla_generate_target_from_python(NAME MTWCodegenThermalPSMCHT_${config}
                    FILE thermalPSMCodegenCHT.py
                    CODEGEN_CFG ${config}
                    OUT_FILES LBMFluidSweep.h LBMFluidSweep.${CODEGEN_FILE_SUFFIX}
                    PSMFluidSweep.h PSMFluidSweep.${CODEGEN_FILE_SUFFIX}
                    LBMConcentrationSweep.h LBMConcentrationSweep.${CODEGEN_FILE_SUFFIX}
                    PSMConcentrationSweep.h PSMConcentrationSweep.${CODEGEN_FILE_SUFFIX}
                    LBMFluidSplitSweep.h LBMFluidSplitSweep.${CODEGEN_FILE_SUFFIX}
                    LBMConcentrationSplitSweep.h LBMConcentrationSplitSweep.${CODEGEN_FILE_SUFFIX}
                    PSMFluidSweepSplit.h PSMFluidSweepSplit.${CODEGEN_FILE_SUFFIX}
                    PSMConcentrationSweepSplit.h PSMConcentrationSweepSplit.${CODEGEN_FILE_SUFFIX}
                    PSMEnergySweep.h PSMEnergySweep.${CODEGEN_FILE_SUFFIX}
                    PackInfoFluid.h PackInfoFluid.${CODEGEN_FILE_SUFFIX}
                    PackInfoConcentration.h PackInfoConcentration.${CODEGEN_FILE_SUFFIX}
                    PackInfoEnergy.h PackInfoEnergy.${CODEGEN_FILE_SUFFIX}
                    InitializeFluidDomain.h InitializeFluidDomain.${CODEGEN_FILE_SUFFIX}
                    InitializeConcentrationDomain.h InitializeConcentrationDomain.${CODEGEN_FILE_SUFFIX}
                    InitializeEnergyDomain.h InitializeEnergyDomain.${CODEGEN_FILE_SUFFIX}
                    BC_Fluid_NoSlip.h BC_Fluid_NoSlip.${CODEGEN_FILE_SUFFIX}
                    BC_Fluid_UBB.h BC_Fluid_UBB.${CODEGEN_FILE_SUFFIX}
                    BC_Fluid_Density.h BC_Fluid_Density.${CODEGEN_FILE_SUFFIX}
                    BC_Fluid_FreeSlip.h BC_Fluid_FreeSlip.${CODEGEN_FILE_SUFFIX}
                    BC_Concentration_NoSlip.h BC_Concentration_NoSlip.${CODEGEN_FILE_SUFFIX}
                    BC_Concentration_UBB.h BC_Concentration_UBB.${CODEGEN_FILE_SUFFIX}
                    BC_Concentration_Density.h BC_Concentration_Density.${CODEGEN_FILE_SUFFIX}
                    BC_Concentration_FreeSlip.h BC_Concentration_FreeSlip.${CODEGEN_FILE_SUFFIX}
                    BC_Concentration_Neumann.h BC_Concentration_Neumann.${CODEGEN_FILE_SUFFIX}
                    BC_Fluid_Outflow.h BC_Fluid_Outflow.${CODEGEN_FILE_SUFFIX}
                    FluidMacroSetter.h FluidMacroSetter.cpp
                    FluidMacroGetter.h FluidMacroGetter.cpp
                    ConcentrationMacroSetter.h ConcentrationMacroSetter.cpp
                    ConcentrationMacroGetter.h ConcentrationMacroGetter.cpp
                    EnergyMacroSetter.h EnergyMacroSetter.cpp
                    EnergyMacroGetter.h EnergyMacroGetter.cpp
                    compute_temperature_field.h compute_temperature_field.cpp
                    GeneralInfoHeader.h
        )
            message(STATUS "generated in thermal PSM is: ${config}" )
            waLBerla_add_executable(NAME thermalPSMCHT_${config} FILES thermalPSMCHT.cpp ../../utilities/InitializerFunctions.cpp
                    DEPENDS walberla::blockforest walberla::core walberla::field walberla::geometry walberla::gpu walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::sqlite walberla::vtk MTWCodegenThermalPSMCHT_${config})
            target_compile_definitions(thermalPSMCHT_${config} PRIVATE Weighting=2)
        endforeach ()


    endif ()
endif ()



