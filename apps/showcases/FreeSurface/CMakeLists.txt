waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME BubblyPoiseuille
      FILES BubblyPoiseuille.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME CapillaryWave
      FILES CapillaryWave.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME DamBreakCylindrical
      FILES DamBreakCylindrical.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME DamBreakRectangular
      FILES DamBreakRectangular.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME DropImpact
      FILES DropImpact.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME DropWetting
      FILES DropWetting.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME GravityWave
      FILES GravityWave.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::postprocessing walberla::timeloop walberla::vtk )

if ( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python( NAME GravityWaveLatticeModelGeneration
         FILE GravityWaveLatticeModelGeneration.py
         OUT_FILES GravityWaveLatticeModel.cpp GravityWaveLatticeModel.h )

   waLBerla_add_executable( NAME GravityWaveCodegen
         FILES GravityWaveCodegen.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm
         walberla::postprocessing walberla::timeloop walberla::vtk GravityWaveLatticeModelGeneration )
endif ()

waLBerla_add_executable( NAME MovingDrop
      FILES MovingDrop.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME RisingBubble
      FILES RisingBubble.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME TaylorBubble
      FILES TaylorBubble.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )