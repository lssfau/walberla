target_sources( walberla_lbm
        PUBLIC
        CellConversionSweep.h
        ConversionFlagsResetSweep.h
        ExcessMassDistributionModel.h
        ExcessMassDistributionSweep.h
        ForceDensitySweep.h
        PdfReconstructionModel.h
        PdfRefillingModel.h
        PdfRefillingSweep.h
        StreamReconstructAdvectSweep.h
        SurfaceDynamicsHandler.h
        ExcessMassDistributionSweep.impl.h
        PdfRefillingSweep.impl.h
)

add_subdirectory( functionality )