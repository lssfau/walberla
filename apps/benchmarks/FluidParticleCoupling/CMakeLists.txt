waLBerla_link_files_to_builddir( "*.dat" )

if( WALBERLA_BUILD_WITH_CODEGEN )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBM FILE GeneratedLBM.py
            OUT_FILES GeneratedLBM.cpp GeneratedLBM.h
            )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBMWithForce FILE GeneratedLBMWithForce.py
            OUT_FILES GeneratedLBMWithForce.cpp GeneratedLBMWithForce.h
            )

   waLBerla_add_executable(NAME SphereWallCollision FILES SphereWallCollision.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling
         walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME SettlingSphereInBox FILES SettlingSphereInBox.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling
         walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME SphereMovingWithPrescribedVelocity FILES SphereMovingWithPrescribedVelocity.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME LubricationForceEvaluation FILES LubricationForceEvaluation.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME DragForceSphere FILES DragForceSphere.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBMWithForce )

   waLBerla_add_executable(NAME ForcesOnSphereNearPlane FILES ForcesOnSphereNearPlane.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME ObliqueWetCollision FILES ObliqueWetCollision.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

   waLBerla_add_executable(NAME MotionSettlingSphere FILES MotionSettlingSphere.cpp
         DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
         walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

else()

waLBerla_add_executable ( NAME SphereWallCollision FILES SphereWallCollision.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME SettlingSphereInBox FILES SettlingSphereInBox.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME SphereMovingWithPrescribedVelocity FILES SphereMovingWithPrescribedVelocity.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME LubricationForceEvaluation FILES LubricationForceEvaluation.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME DragForceSphere FILES DragForceSphere.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME ForcesOnSphereNearPlane FILES ForcesOnSphereNearPlane.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME ObliqueWetCollision FILES ObliqueWetCollision.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

endif()

waLBerla_add_executable ( NAME ObliqueDryCollision FILES ObliqueDryCollision.cpp
      DEPENDS walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing )
