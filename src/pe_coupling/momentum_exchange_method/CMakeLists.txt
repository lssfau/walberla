target_sources( walberla_pe_coupling
    PRIVATE
    BodyMapping.h     
    )
add_subdirectory( boundary )
add_subdirectory( destruction )
add_subdirectory( restoration )
