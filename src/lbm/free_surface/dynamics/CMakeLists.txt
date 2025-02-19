target_sources( walberla_lbm
        PRIVATE
        CellConversionSweep.h
        ConversionFlagsResetSweep.h
        ExcessMassDistributionModel.h
        ExcessMassDistributionSweep.h
        ExcessMassDistributionSweep.impl.h
        ForceDensitySweep.h
        PdfReconstructionModel.h
        PdfRefillingModel.h
        PdfRefillingSweep.h
        PdfRefillingSweep.impl.h
        StreamReconstructAdvectSweep.h
        SurfaceDynamicsHandler.h
        )

add_subdirectory( functionality )