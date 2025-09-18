target_sources( walberla_geometry
        PUBLIC
        RGBAImage.h
        BasicVoxelFileReader.impl.h
        GrayScaleImage.h
        BasicVoxelFileReader.h
        BinaryRawFile.h
        VoxelFileReader.impl.h
        VoxelFileReader.h
)

target_sources( walberla_geometry
        PRIVATE
        RGBAImage.cpp
        GrayScaleImage.cpp
        BinaryRawFile.cpp
)
