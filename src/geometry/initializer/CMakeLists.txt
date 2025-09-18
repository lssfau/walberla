target_sources( walberla_geometry
        PUBLIC
        BoundaryFromCellInterval.h
        BoundaryFromBody.impl.h
        BoundaryFromVoxelFile.impl.h
        BoundaryFromBody.h
        ScalarFieldFromBody.h
        BoundaryFromVoxelFile.h
        ScalarFieldFromBody.impl.h
        ScalarFieldFromGrayScaleImage.h
        ScalarFieldFromCellInterval.h
        BoundaryFromImage.h
        InitializationManager.h
        BoundaryFromDomainBorder.impl.h
        BoundaryFromDomainBorder.h
        BoundaryFromImage.impl.h
        Initializer.h
        OverlapFieldFromBody.h
        BoundarySetter.h
        ScalarFieldFromCellInterval.impl.h
        BoundaryFromCellInterval.impl.h
        BoundarySetterFlagFieldSpecialization.h
)

target_sources( walberla_geometry
        PRIVATE
        BoundaryFromVoxelFile.cpp
        ScalarFieldFromGrayScaleImage.cpp
        InitializationManager.cpp
        OverlapFieldFromBody.cpp
)
