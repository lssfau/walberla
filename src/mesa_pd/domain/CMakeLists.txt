target_sources( walberla_mesa_pd
        PUBLIC
        IDomain.h
        InfiniteDomain.h
        BlockForestDataHandling.h
        InfoCollection.h
        BlockForestDomain.h
)

target_sources( walberla_mesa_pd
        PRIVATE
        BlockForestDomain.cpp
        BlockForestDataHandling.cpp
)
