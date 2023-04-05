waLBerla_link_files_to_builddir( *.prm )

if( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python( NAME      AntidunesLatticeModelGeneration
                                         FILE      AntidunesLatticeModelGeneration.py
                                         OUT_FILES AntidunesLatticeModel.cpp AntidunesLatticeModel.h )

   waLBerla_add_executable(NAME     Antidunes
                           FILES    Antidunes.cpp PIDController.cpp
                           DEPENDS  blockforest boundary core domain_decomposition field lbm mesa_pd lbm_mesapd_coupling
                                    postprocessing timeloop vtk AntidunesLatticeModelGeneration)
endif()

waLBerla_add_executable(NAME     BedGeneration
        FILES    BedGeneration.cpp
        DEPENDS  blockforest core domain_decomposition mesa_pd vtk)
