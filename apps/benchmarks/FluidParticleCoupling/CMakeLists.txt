waLBerla_link_files_to_builddir( "*.dat" )

if( WALBERLA_BUILD_WITH_CODEGEN )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBM FILE GeneratedLBM.py
            OUT_FILES GeneratedLBM.cpp GeneratedLBM.h
            )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBMWithForce FILE GeneratedLBMWithForce.py
            OUT_FILES GeneratedLBMWithForce.cpp GeneratedLBMWithForce.h
            )

   waLBerla_add_executable(NAME SphereWallCollision FILES SphereWallCollision.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling
         mesa_pd postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

   waLBerla_add_executable(NAME SettlingSphereInBox FILES SettlingSphereInBox.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling
         mesa_pd postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

   waLBerla_add_executable(NAME SphereMovingWithPrescribedVelocity FILES SphereMovingWithPrescribedVelocity.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm mesa_pd lbm_mesapd_coupling
         postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

   waLBerla_add_executable(NAME LubricationForceEvaluation FILES LubricationForceEvaluation.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd
         postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

   waLBerla_add_executable(NAME DragForceSphere FILES DragForceSphere.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd
         postprocessing timeloop vtk FluidParticleCouplingGeneratedLBMWithForce)

   waLBerla_add_executable(NAME ForcesOnSphereNearPlane FILES ForcesOnSphereNearPlane.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd
         postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

   waLBerla_add_executable(NAME ObliqueWetCollision FILES ObliqueWetCollision.cpp
         DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd
         postprocessing timeloop vtk FluidParticleCouplingGeneratedLBM)

else()

waLBerla_add_executable ( NAME SphereWallCollision FILES SphereWallCollision.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME SettlingSphereInBox FILES SettlingSphereInBox.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME SphereMovingWithPrescribedVelocity FILES SphereMovingWithPrescribedVelocity.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm mesa_pd lbm_mesapd_coupling postprocessing timeloop vtk )

waLBerla_add_executable ( NAME LubricationForceEvaluation FILES LubricationForceEvaluation.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME DragForceSphere FILES DragForceSphere.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME ForcesOnSphereNearPlane FILES ForcesOnSphereNearPlane.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME ObliqueWetCollision FILES ObliqueWetCollision.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

endif()



waLBerla_add_executable ( NAME ObliqueDryCollision FILES ObliqueDryCollision.cpp
      DEPENDS blockforest core mesa_pd postprocessing )