target_sources( blockforest
      PRIVATE
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
      Cartesian.cpp
      DynamicParMetis.cpp
      StaticCurve.cpp
      StaticParMetis.cpp
      )


add_subdirectory( weight_assignment )
add_subdirectory( level_determination )