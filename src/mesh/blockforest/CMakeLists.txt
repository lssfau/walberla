target_sources( walberla_mesh
        PUBLIC
        RefinementSelection.h
        BlockWorkloadMemory.h
        BlockExclusion.h
        BlockForestInitialization.h
)

target_sources( walberla_mesh
        PRIVATE
        BlockForestInitialization.cpp
)