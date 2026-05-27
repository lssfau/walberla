if(WALBERLA_BUILD_WITH_METIS)
    if(METIS_FOUND)
        target_link_libraries (walberla_core PUBLIC ${METIS_LIBRARIES})
    endif()
    target_sources (
        walberla_core
        PRIVATE MetisWrapper.cpp
        PUBLIC MetisWrapper.h
    )
endif()

if(WALBERLA_BUILD_WITH_PARMETIS)
    if(NOT "${PARMETIS_LIBRARY}" STREQUAL "")
        target_link_libraries (walberla_core PUBLIC ${PARMETIS_LIBRARY})
    endif()
    target_sources (
        walberla_core
        PRIVATE ParMetisWrapper.cpp
        PUBLIC ParMetisWrapper.h
    )
endif()
