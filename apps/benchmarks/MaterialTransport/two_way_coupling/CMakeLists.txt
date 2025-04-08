####################################################################################
## to run this in 2 dimensions, set the flag -Duse_2D = ON in the CMake settings ###
## to run this in 3 dimensions, set the flag -Duse_2D = OFF in the CMake settings ##
####################################################################################

waLBerla_link_files_to_builddir("*.prm")

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (use_2D)
        set(dimensions 2)
        add_definitions(-Duse_2D)
    else ()
        set(dimensions 3)
    endif()
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))

        set(config "dimension_${dimensions}")
        message(STATUS "Config value is: ${config}" )
        waLBerla_generate_target_from_python(NAME MTW_2WayCodegen
                FILE MTW_2WayCodegen.py
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
                #BC_Concentration_UBB.h BC_Concentration_UBB.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_Density.h BC_Concentration_Density.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_FreeSlip.h BC_Concentration_FreeSlip.${CODEGEN_FILE_SUFFIX}
                BC_Concentration_Neumann.h BC_Concentration_Neumann.${CODEGEN_FILE_SUFFIX}
                FluidMacroSetter.h FluidMacroSetter.cpp
                FluidMacroGetter.h FluidMacroGetter.cpp
                ConcentrationMacroSetter.h ConcentrationMacroSetter.cpp
                ConcentrationMacroGetter.h ConcentrationMacroGetter.cpp
                GeneralInfoHeader.h
        )

        waLBerla_add_executable(NAME MTW_2Way FILES MaterialTransport.cpp ../utilities/InitializerFunctions.cpp
                DEPENDS walberla::blockforest walberla::core walberla::field walberla::geometry walberla::gpu walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::sqlite walberla::vtk MTW_2WayCodegen)
        target_compile_definitions(MTW_2Way PRIVATE Weighting=2)
    endif ()
endif ()
