if(WALBERLA_BUILD_WITH_METIS)
    target_link_libraries (walberla_core PUBLIC ${METIS_LIBRARIES})

    target_sources (
        walberla_core
        PRIVATE MetisWrapper.cpp
        PUBLIC MetisWrapper.h
    )
endif()

if(WALBERLA_BUILD_WITH_PARMETIS)
    target_link_libraries (walberla_core PUBLIC ${PARMETIS_LIBRARY})

    target_sources (
        walberla_core
        PRIVATE ParMetisWrapper.cpp
        PUBLIC ParMetisWrapper.h
    )
endif()
