target_sources( lbm
        PRIVATE
        CellConversionSweep.h
        ConversionFlagsResetSweep.h
        ExcessMassDistributionModel.h
        ExcessMassDistributionSweep.h
        ExcessMassDistributionSweep.impl.h
        ForceWeightingSweep.h
        PdfReconstructionModel.h
        PdfRefillingModel.h
        PdfRefillingSweep.h
        PdfRefillingSweep.impl.h
        StreamReconstructAdvectSweep.h
        SurfaceDynamicsHandler.h
        )

add_subdirectory( functionality )