target_sources( walberla_geometry
        PUBLIC
        Torus.h
        BodyFromConfig.h
        Ellipsoid.h
        Sphere.h
        BodyLogic.h
        AABBBody.h
        Cylinder.h
        BodyOverlapFunctions.h
        BodyOverlapFunctions.impl.h
        DynamicBody.h
)

target_sources( walberla_geometry
        PRIVATE
        BodyFromConfig.cpp
        Cylinder.cpp
        Ellipsoid.cpp
        Torus.cpp
        Sphere.cpp
)
