waLBerla_link_files_to_builddir( *.prm )

if ( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python( NAME AntidunesLatticeModelGeneration
         FILE AntidunesLatticeModelGeneration.py
         OUT_FILES AntidunesLatticeModel.cpp AntidunesLatticeModel.h )

   waLBerla_add_executable( NAME Antidunes
         FILES Antidunes.cpp PIDController.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field
         walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling walberla::postprocessing walberla::timeloop
         walberla::vtk AntidunesLatticeModelGeneration )
endif ()

waLBerla_add_executable( NAME BedGeneration
      FILES BedGeneration.cpp
      DEPENDS walberla::blockforest walberla::core walberla::domain_decomposition walberla::mesa_pd walberla::vtk )
