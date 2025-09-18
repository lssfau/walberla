target_sources( walberla_lbm
        PUBLIC
        ContactAngle.h
        CurvatureModel.h
        CurvatureModel.impl.h
        CurvatureSweep.h
        CurvatureSweep.impl.h
        DetectWettingSweep.h
        ExtrapolateNormalsSweep.h
        ExtrapolateNormalsSweep.impl.h
        NormalSweep.h
        NormalSweep.impl.h
        ObstacleFillLevelSweep.h
        ObstacleFillLevelSweep.impl.h
        ObstacleNormalSweep.h
        ObstacleNormalSweep.impl.h
        SmoothingSweep.h
        SmoothingSweep.impl.h
        SurfaceGeometryHandler.h
        Utility.h
)

target_sources( walberla_lbm
        PRIVATE
        Utility.cpp
)
