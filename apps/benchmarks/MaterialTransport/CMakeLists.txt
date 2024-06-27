waLBerla_link_files_to_builddir("*.prm")

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))

        waLBerla_generate_target_from_python(NAME MTWCodegen
                FILE MTWCodegen.py
                OUT_FILES LBMFluidSweep.h LBMFluidSweep.cpp
                LBMConcentrationSweep.h LBMConcentrationSweep.cpp
                LBMFluidSplitSweep.h LBMFluidSplitSweep.cpp
                LBMConcentrationSplitSweep.h LBMConcentrationSplitSweep.cpp
                PackInfoFluid.h PackInfoFluid.cpp
                PackInfoConcentration.h PackInfoConcentration.cpp
                InitializeFluidDomain.h InitializeFluidDomain.cpp
                InitializeConcentrationDomain.h InitializeConcentrationDomain.cpp
                BC_Fluid_NoSlip.h BC_Fluid_NoSlip.cpp
                BC_Fluid_UBB.h BC_Fluid_UBB.cpp
                BC_Fluid_Density.h BC_Fluid_Density.cpp
                BC_Fluid_FreeSlip.h BC_Fluid_FreeSlip.cpp
                BC_Concentration_NoSlip.h BC_Concentration_NoSlip.cpp
                BC_Concentration_UBB.h BC_Concentration_UBB.cpp
                BC_Concentration_Density.h BC_Concentration_Density.cpp
                BC_Concentration_FreeSlip.h BC_Concentration_FreeSlip.cpp
                FluidMacroSetter.h FluidMacroSetter.cpp
                FluidMacroGetter.h FluidMacroGetter.cpp
                ConcentrationMacroSetter.h ConcentrationMacroSetter.cpp
                ConcentrationMacroGetter.h ConcentrationMacroGetter.cpp
                GeneralInfoHeader.h
        )

        waLBerla_add_executable(NAME MaterialTransport FILES MaterialTransport.cpp
                DEPENDS blockforest core field geometry gpu lbm lbm_mesapd_coupling mesa_pd sqlite vtk MTWCodegen)
        target_compile_definitions(MaterialTransport PRIVATE Weighting=2)
    endif ()
endif ()
