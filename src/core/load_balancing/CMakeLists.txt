if(WALBERLA_BUILD_WITH_METIS)
    target_link_libraries (walberla_core PUBLIC ${METIS_LIBRARIES})
endif()

if(WALBERLA_BUILD_WITH_PARMETIS)
    target_link_libraries (walberla_core PUBLIC ${PARMETIS_LIBRARY})
endif()

target_sources (
    walberla_core
    PRIVATE MetisWrapper.cpp ParMetisWrapper.cpp
    PUBLIC MetisWrapper.h ParMetisWrapper.h
)
