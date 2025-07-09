add_subdirectory( FluidizedBed )
add_subdirectory( FreeSurface )
add_subdirectory( HeatConduction )
add_subdirectory( LightRisingParticleInFluidAMR )
add_subdirectory( Mixer )
add_subdirectory( ParticlePacking )
add_subdirectory( PegIntoSphereBed )
add_subdirectory( Piping )
add_subdirectory( ThermalFluidizedBed )
if ( WALBERLA_BUILD_WITH_CODEGEN)

   add_subdirectory( Antidunes )
   add_subdirectory( FlowAroundSphere )
   add_subdirectory( KelvinHelmholtzInstability )
   
   if (WALBERLA_BUILD_WITH_PYTHON)
      add_subdirectory( PhaseFieldAllenCahn )
   endif ()

   if (WALBERLA_DOUBLE_ACCURACY)
      add_subdirectory( Thermocapillary )
   endif ()

endif()

