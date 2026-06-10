add_subdirectory (CNT)
add_subdirectory (GranularGas)
add_subdirectory (IntegratorAccuracy)
add_subdirectory (MeshDistance)
add_subdirectory (ProbeVsExtraMessage)

if(WALBERLA_BUILD_LEGACY_LBM)
    add_subdirectory (ComplexGeometry)
    add_subdirectory (CouetteFlow)
    add_subdirectory (FluidParticleCoupling)
    add_subdirectory (FluidParticleCouplingWithLoadBalancing)
    add_subdirectory (FreeSurfaceAdvection)
    add_subdirectory (NonUniformGrid)
    add_subdirectory (PoiseuilleChannel)
    add_subdirectory (SchaeferTurek)
    add_subdirectory (UniformGrid)
    if(WALBERLA_BUILD_WITH_DOUBLE_ACCURACY)
        # The vectorized PSM codes are only correctly generated for double precision.
        add_subdirectory (FluidizedBed)
        add_subdirectory (Percolation)
    endif()
endif()

if(WALBERLA_BUILD_WITH_PYTHON)
    add_subdirectory (FieldCommunication)

    if(WALBERLA_BUILD_WITH_CODEGEN AND WALBERLA_BUILD_LEGACY_LBM)
        add_subdirectory (NonUniformGridCPU)
        add_subdirectory (PhaseFieldAllenCahn)
        add_subdirectory (TurbulentChannel)
        add_subdirectory (UniformGridCPU)
    endif()

    if(WALBERLA_BUILD_WITH_CODEGEN
       AND WALBERLA_BUILD_WITH_GPU_SUPPORT
       AND WALBERLA_BUILD_LEGACY_LBM
    )
        add_subdirectory (NonUniformGridGPU)
        add_subdirectory (UniformGridGPU)
    endif()
endif()
