waLBerla_link_files_to_builddir( "*.dat" )

if( WALBERLA_BUILD_WITH_CODEGEN )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBM FILE GeneratedLBM.py
            OUT_FILES GeneratedLBM.cpp GeneratedLBM.h
            )

    waLBerla_generate_target_from_python(NAME FluidParticleCouplingGeneratedLBMWithForce FILE GeneratedLBMWithForce.py
            OUT_FILES GeneratedLBMWithForce.cpp GeneratedLBMWithForce.h
            )

    add_executable ( SphereWallCollision SphereWallCollision.cpp )
    target_link_libraries ( SphereWallCollision
            walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field
            walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing walberla::timeloop
            walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable ( SettlingSphereInBox SettlingSphereInBox.cpp )
    target_link_libraries ( SettlingSphereInBox walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable ( SphereMovingWithPrescribedVelocity SphereMovingWithPrescribedVelocity.cpp )
    target_link_libraries ( SphereMovingWithPrescribedVelocity walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable ( LubricationForceEvaluation LubricationForceEvaluation.cpp )
    target_link_libraries ( LubricationForceEvaluation walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable( DragForceSphere DragForceSphere.cpp )
    target_link_libraries ( DragForceSphere walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBMWithForce )

    add_executable ( ForcesOnSphereNearPlane ForcesOnSphereNearPlane.cpp )
    target_link_libraries ( ForcesOnSphereNearPlane walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable ( ObliqueWetCollision ObliqueWetCollision.cpp )
    target_link_libraries ( ObliqueWetCollision walberla::blockforest walberla::boundary walberla::core
            walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
            walberla::postprocessing walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

    add_executable( MotionSettlingSphere MotionSettlingSphere.cpp )
    target_link_libraries ( MotionSettlingSphere walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
            walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing
            walberla::timeloop walberla::vtk FluidParticleCouplingGeneratedLBM )

else()

add_executable ( SphereWallCollision SphereWallCollision.cpp )
target_link_libraries ( SphereWallCollision walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing
        walberla::timeloop walberla::vtk )

add_executable ( SettlingSphereInBox SettlingSphereInBox.cpp )
target_link_libraries ( SettlingSphereInBox walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing
        walberla::timeloop walberla::vtk )

add_executable ( SphereMovingWithPrescribedVelocity SphereMovingWithPrescribedVelocity.cpp )
target_link_libraries ( SphereMovingWithPrescribedVelocity walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling walberla::postprocessing
        walberla::timeloop walberla::vtk )

add_executable ( LubricationForceEvaluation LubricationForceEvaluation.cpp )
target_link_libraries ( LubricationForceEvaluation walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
        walberla::postprocessing walberla::timeloop walberla::vtk )

add_executable ( DragForceSphere DragForceSphere.cpp )
target_link_libraries ( DragForceSphere walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling
        walberla::mesa_pd walberla::postprocessing walberla::timeloop walberla::vtk )

add_executable ( ForcesOnSphereNearPlane ForcesOnSphereNearPlane.cpp )
target_link_libraries ( ForcesOnSphereNearPlane walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
        walberla::postprocessing walberla::timeloop walberla::vtk )

add_executable ( ObliqueWetCollision ObliqueWetCollision.cpp )
target_link_libraries ( ObliqueWetCollision walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
        walberla::postprocessing walberla::timeloop walberla::vtk )

endif()

add_executable ( ObliqueDryCollision ObliqueDryCollision.cpp )
target_link_libraries ( ObliqueDryCollision walberla::blockforest walberla::core walberla::mesa_pd
        walberla::postprocessing )
