waLBerla_link_files_to_builddir( *.prm )

add_executable ( BubblyPoiseuille BubblyPoiseuille.cpp )
target_link_libraries ( BubblyPoiseuille walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( CapillaryWave CapillaryWave.cpp )
target_link_libraries ( CapillaryWave walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( DamBreakCylindrical DamBreakCylindrical.cpp )
target_link_libraries ( DamBreakCylindrical walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( DamBreakRectangular DamBreakRectangular.cpp )
target_link_libraries ( DamBreakRectangular walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( DropImpact DropImpact.cpp )
target_link_libraries ( DropImpact walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( DropWetting DropWetting.cpp )
target_link_libraries ( DropWetting walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( GravityWave GravityWave.cpp )
target_link_libraries ( GravityWave walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

if ( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python( NAME GravityWaveLatticeModelGeneration
         FILE GravityWaveLatticeModelGeneration.py
         OUT_FILES GravityWaveLatticeModel.cpp GravityWaveLatticeModel.h )

   add_executable ( GravityWaveCodegen GravityWaveCodegen.cpp )
   target_link_libraries ( GravityWaveCodegen walberla::blockforest walberla::boundary walberla::core
           walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
           walberla::vtk GravityWaveLatticeModelGeneration )
endif ()

add_executable ( MovingDrop MovingDrop.cpp )
target_link_libraries ( MovingDrop walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( RisingBubble RisingBubble.cpp )
target_link_libraries ( RisingBubble walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( TaylorBubble TaylorBubble.cpp )
target_link_libraries ( TaylorBubble walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )