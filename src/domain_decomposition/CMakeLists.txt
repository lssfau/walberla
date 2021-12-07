###################################################################################################
#
# Module domain_decomposition
#
###################################################################################################

add_library( domain_decomposition )
target_link_libraries( domain_decomposition PUBLIC core )
target_sources( domain_decomposition
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