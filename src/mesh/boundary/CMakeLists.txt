target_sources( walberla_mesh
        PUBLIC
        ColorToBoundaryMapper.h
        BoundaryLocationFunction.h
        BoundaryLocation.h
        BoundaryUIDFaceDataSource.h
        BoundarySetup.h
        BoundaryInfo.h
)

target_sources( walberla_mesh
        PRIVATE
        BoundarySetup.cpp
)
