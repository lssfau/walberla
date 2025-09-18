OPTION( WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID "BlockForest: enables/disables BlockID implementation as a primitive data type (uint_t)" ON )

mark_as_advanced( WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID )

configure_file( CMakeDefs.in.h CMakeDefs.h )

add_library( walberla_blockforest )
add_library( walberla::blockforest ALIAS walberla_blockforest )
target_link_libraries( walberla_blockforest PUBLIC walberla::communication walberla::core walberla::domain_decomposition walberla::stencil )
target_sources( walberla_blockforest
        PUBLIC
        AABBRefinementSelection.h
        all.h
        BlockDataHandling.h
        BlockForestEvaluation.h
        BlockForestFile.h
        BlockForest.h
        Block.h
        BlockID.h
        BlockNeighborhoodConstruction.h
        BlockNeighborhoodSection.h
        BlockReconstruction.h
        CMakeDefs.in.h
        GlobalLoadBalancing.h
        HilbertCurveConstruction.h
        Initialization.h
        OutputColor.h
        PhantomBlockForest.h
        PhantomBlock.h
        SetupBlockForest.h
        SetupBlock.h
        StructuredBlockForest.h
        Types.h
        Utility.h
)

target_sources( walberla_blockforest
        PRIVATE
        Block.cpp
        BlockDataHandling.cpp
        BlockForest.cpp
        BlockForestEvaluation.cpp
        BlockID.cpp
        BlockNeighborhoodConstruction.cpp
        BlockReconstruction.cpp
        Initialization.cpp
        OutputColor.cpp
        PhantomBlock.cpp
        PhantomBlockForest.cpp
        SetupBlock.cpp
        SetupBlockForest.cpp
        StructuredBlockForest.cpp
)

add_subdirectory( communication )
add_subdirectory( loadbalancing )
