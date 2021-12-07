target_sources( lbm_mesapd_coupling
    PRIVATE
    SubCyclingManager.h
    AddHydrodynamicInteractionKernel.h
    AddAccelerationOnParticlesKernel.h
    SubCyclingManager.cpp
    ParticleFunctions.h
    ResetHydrodynamicForceTorqueKernel.h
    ParticleSelector.h
    LubricationCorrectionKernel.h
    AddForceOnParticlesKernel.h
    OmegaBulkAdaption.h
    AverageHydrodynamicForceTorqueKernel.h
    InspectionProbe.h
    InitializeHydrodynamicForceTorqueForAveragingKernel.h     
    )
