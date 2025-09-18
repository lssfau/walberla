target_sources( walberla_mesa_pd
        PUBLIC
        AnalyticCollisionFunctions.h
        GJK.h
        GeneralContactDetection.h
        Support.h
        EPA.h
        BroadPhase.h
        AnalyticContactDetection.h
)

target_sources( walberla_mesa_pd
        PRIVATE
        GJK.cpp
        EPA.cpp
)
