target_sources( walberla_pe_coupling
    PRIVATE
    InfoCollection.cpp
    BlockInfo.h
    InfoCollection.h     
    )
add_subdirectory( level_determination )
add_subdirectory( weight_assignment )
