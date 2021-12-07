target_sources( lbm_mesapd_coupling
    PRIVATE
    BlockInfo.h
    InfoCollection.h     
    )

add_subdirectory( level_determination )
add_subdirectory( weight_assignment )