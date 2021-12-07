target_sources( lbm
    PRIVATE
    MomentumDensity.h
    MacroscopicValueCalculation.h
    Equilibrium.h
    DensityAndVelocity.h
    QCriterion.h
    DensityVelocityCallback.h
    AddToStorage.h
    QCriterionFieldWriter.h
    VelocityFieldWriter.h
    ShearRate.h
    PdfField.h
    Density.h
    Vorticity.h
    DensityAndMomentumDensity.h
    Adaptors.h
    PressureTensor.h     
    )
add_subdirectory( initializer )