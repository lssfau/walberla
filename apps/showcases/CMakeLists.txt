if( WALBERLA_BUILD_LEGACY_LBM )
   add_subdirectory( FluidizedBed )
   add_subdirectory( FreeSurface )
endif()
add_subdirectory( HeatConduction )
if( WALBERLA_BUILD_LEGACY_LBM )
   add_subdirectory( LightRisingParticleInFluidAMR )
endif()
add_subdirectory( Mixer )
add_subdirectory( ParticlePacking )
add_subdirectory( PegIntoSphereBed )
if( WALBERLA_BUILD_LEGACY_LBM )
   add_subdirectory( Piping )
endif()
if ( WALBERLA_BUILD_WITH_CODEGEN)

   if( WALBERLA_BUILD_LEGACY_LBM )
      add_subdirectory( Antidunes )
      add_subdirectory( FlowAroundSphere )
   endif()
   
   if (WALBERLA_BUILD_WITH_PYTHON AND WALBERLA_BUILD_LEGACY_LBM)
      add_subdirectory( PhaseFieldAllenCahn )

      add_subdirectory( KelvinHelmholtzInstability )

      if (WALBERLA_DOUBLE_ACCURACY)
         add_subdirectory( Thermocapillary )
      endif ()
   endif ()

endif()

