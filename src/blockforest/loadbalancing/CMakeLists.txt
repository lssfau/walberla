target_sources( walberla_blockforest
      PUBLIC
      BlockInfo.h
      Cartesian.h
      DynamicCurve.h
      DynamicDiffusive.h
      DynamicParMetis.h
      InfoCollection.h
      NoPhantomData.h
      PODPhantomData.h
      StaticCurve.h
      StaticParMetis.h
)

target_sources( walberla_blockforest
        PRIVATE
        Cartesian.cpp
        DynamicParMetis.cpp
        StaticCurve.cpp
        StaticParMetis.cpp
)


add_subdirectory( weight_assignment )
add_subdirectory( level_determination )