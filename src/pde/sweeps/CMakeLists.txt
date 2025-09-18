target_sources( walberla_pde
        PUBLIC
        StencilSweepBase.h
        SweepBase.h
        RBGS.h
        RBGSFixedStencil.h
        SOR.h
        SORFixedStencil.h
        Jacobi.h
        StencilFieldSweepBase.h
        Multigrid.h
        Multigrid.impl.h
        JacobiFixedStencil.h
)

target_sources( walberla_pde
        PRIVATE
        SweepBase.cpp
)
