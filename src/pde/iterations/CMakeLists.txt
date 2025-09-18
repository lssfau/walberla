target_sources( walberla_pde
        PUBLIC
        JacobiIteration.h
        RBGSIteration.h
        CGFixedStencilIteration.h
        VCycles.h
        VCycles.impl.h
        CGIteration.h
)

target_sources( walberla_pde
        PRIVATE
        JacobiIteration.cpp
        RBGSIteration.cpp
)
