target_sources( pe
    PRIVATE
    InfoCollection.cpp
    InfoCollection.h     
    )
add_subdirectory( level_determination )
add_subdirectory( weight_assignment )