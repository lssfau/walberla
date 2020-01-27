waLBerla_link_files_to_builddir( "*.dat" )

waLBerla_python_file_generates(GeneratedLBM.py
      GeneratedLBM.cpp GeneratedLBM.h
      )

waLBerla_python_file_generates(GeneratedLBMWithForce.py
      GeneratedLBMWithForce.cpp GeneratedLBMWithForce.h
      )

waLBerla_add_executable ( NAME SphereWallCollision FILES SphereWallCollision.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME SettlingSphereTenCate FILES SettlingSphereInBox.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME SphereMovingWithPrescribedVelocity FILES SphereMovingWithPrescribedVelocity.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm mesa_pd lbm_mesapd_coupling postprocessing timeloop vtk )

waLBerla_add_executable ( NAME LubricationForceEvaluation FILES LubricationForceEvaluation.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME DragForceSphere FILES DragForceSphere.cpp GeneratedLBMWithForce.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME ForcesOnSphereNearPlane FILES ForcesOnSphereNearPlane.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )

waLBerla_add_executable ( NAME ObliqueDryCollision FILES ObliqueDryCollision.cpp
                          DEPENDS blockforest core mesa_pd postprocessing )

waLBerla_add_executable ( NAME ObliqueWetCollision FILES ObliqueWetCollision.cpp GeneratedLBM.py
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )