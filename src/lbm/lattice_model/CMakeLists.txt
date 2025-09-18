target_sources( walberla_lbm
        PUBLIC
        D3Q27.h
        CollisionModel.h
        LatticeModelBase.h
        ForceModel.h
        D3Q19.h
        EquilibriumDistribution.h
        D3Q15.h
        NeighborsStencil.h
        D2Q9.h
        SmagorinskyLES.h
)

target_sources( walberla_lbm
        PRIVATE
        CollisionModel.cpp
)
