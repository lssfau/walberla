target_sources( pe_coupling
    PRIVATE
    BodyMapping.h     
    )
add_subdirectory( boundary )
add_subdirectory( destruction )
add_subdirectory( restoration )
