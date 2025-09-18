if ( METIS_FOUND )
    target_link_libraries( walberla_core PUBLIC ${METIS_LIBRARIES} )
endif()

if ( NOT "${PARMETIS_LIBRARY}" STREQUAL "" )
    target_link_libraries( walberla_core PUBLIC ${PARMETIS_LIBRARY} )
endif()

target_sources( walberla_core
        PRIVATE
        MetisWrapper.cpp
        ParMetisWrapper.cpp
)

target_sources( walberla_core
        PUBLIC
        MetisWrapper.h
        ParMetisWrapper.h
)