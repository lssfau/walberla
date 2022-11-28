waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable(NAME    BubblyPoiseuille
                        FILES   BubblyPoiseuille.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    CapillaryWave
                        FILES   CapillaryWave.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    DamBreakCylindrical
                        FILES   DamBreakCylindrical.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    DamBreakRectangular
                        FILES   DamBreakRectangular.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    DropImpact
                        FILES   DropImpact.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    DropWetting
                        FILES   DropWetting.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    GravityWave
                        FILES   GravityWave.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

if( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python( NAME      GravityWaveLatticeModelGeneration
                                         FILE      GravityWaveLatticeModelGeneration.py
                                         OUT_FILES GravityWaveLatticeModel.cpp GravityWaveLatticeModel.h )

   waLBerla_add_executable(NAME    GravityWaveCodegen
                           FILES   GravityWaveCodegen.cpp
                           DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk
                                   GravityWaveLatticeModelGeneration)
endif()

waLBerla_add_executable(NAME    RisingBubble
                        FILES   RisingBubble.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    TaylorBubble
                        FILES   TaylorBubble.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    ZalesakDisk
                        FILES   ZalesakDisk.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)