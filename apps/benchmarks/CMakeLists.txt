add_subdirectory( AdaptiveMeshRefinementFluidParticleCoupling )
add_subdirectory( ComplexGeometry )
add_subdirectory( DEM )
add_subdirectory( MeshDistance )
add_subdirectory( CouetteFlow )
add_subdirectory( FluidParticleCoupling )
add_subdirectory( FluidParticleCouplingWithLoadBalancing )
add_subdirectory( ForcesOnSphereNearPlaneInShearFlow )
add_subdirectory( GranularGas )
add_subdirectory( IntegratorAccuracy )
add_subdirectory( LennardJones )
add_subdirectory( NonUniformGrid )
add_subdirectory( MotionSingleHeavySphere )
add_subdirectory( PoiseuilleChannel )
add_subdirectory( ProbeVsExtraMessage )
add_subdirectory( SchaeferTurek )
add_subdirectory( UniformGrid )

if ( WALBERLA_BUILD_WITH_PYTHON )
   add_subdirectory( FieldCommunication )

   if ( WALBERLA_BUILD_WITH_CODEGEN )
      add_subdirectory( UniformGridGenerated )
      add_subdirectory( PhaseFieldAllenCahn )
   endif()

   if ( WALBERLA_BUILD_WITH_CODEGEN AND WALBERLA_BUILD_WITH_CUDA )
      add_subdirectory( UniformGridGPU )
   endif()

endif()



