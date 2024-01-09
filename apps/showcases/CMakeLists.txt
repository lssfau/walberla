add_subdirectory( BidisperseFluidizedBed )
add_subdirectory( CombinedResolvedUnresolved )
add_subdirectory( FluidizedBed )
add_subdirectory( FreeSurface )
add_subdirectory( HeatConduction )
add_subdirectory( LightRisingParticleInFluidAMR )
add_subdirectory( Mixer )
add_subdirectory( ParticlePacking )
add_subdirectory( PegIntoSphereBed )
if ( WALBERLA_BUILD_WITH_CODEGEN)

   add_subdirectory( Antidunes )

   if (WALBERLA_BUILD_WITH_PYTHON)
      add_subdirectory( PhaseFieldAllenCahn )
   endif ()

   if(NOT WALBERLA_BUILD_WITH_OPENMP)
      add_subdirectory( PorousMedia )
   endif()

   if (WALBERLA_DOUBLE_ACCURACY)
      add_subdirectory( Thermocapillary )
   endif ()

endif()

