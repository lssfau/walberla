target_sources( lbm
        PRIVATE
        BlockStateDetectorSweep.h
        FlagDefinitions.h
        FlagInfo.h
        FlagInfo.impl.h
        InitFunctions.h
        InterfaceFromFillLevel.h
        LoadBalancing.h
        MaxVelocityComputer.h
        SurfaceMeshWriter.h
        TotalMassComputer.h
        VtkWriter.h
        )

add_subdirectory( boundary )
add_subdirectory( bubble_model )
add_subdirectory( dynamics )
add_subdirectory( surface_geometry )
