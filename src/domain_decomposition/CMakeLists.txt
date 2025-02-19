###################################################################################################
#
# Module domain_decomposition
#
###################################################################################################

add_library( walberla_domain_decomposition )
add_library( walberla::domain_decomposition ALIAS walberla_domain_decomposition )
target_link_libraries( walberla_domain_decomposition PUBLIC walberla::core )
target_sources( walberla_domain_decomposition
      PRIVATE
      BlockCounter.h
      BlockDataHandling.h
      BlockDataID.h
      BlockStorage.h
      BlockSweepWrapper.h
      IBlock.h
      IBlockID.h
      MakeBlockDataInitFunction.h
      MapPointToPeriodicDomain.h
      PeriodicIntersect.h
      PeriodicIntersectionVolume.h
      SharedSweep.h
      StructuredBlockStorageCellMapping.h
      StructuredBlockStorage.h
      BlockDataHandling.cpp
      BlockStorage.cpp
      IBlock.cpp
      MapPointToPeriodicDomain.cpp
      PeriodicIntersect.cpp
      PeriodicIntersectionVolume.cpp
      StructuredBlockStorage.cpp
      )

###################################################################################################                        