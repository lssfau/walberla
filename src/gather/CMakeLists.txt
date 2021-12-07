###################################################################################################
#
# Module gather ( all-to-one communication for output )
#
###################################################################################################

add_library( gather )
target_link_libraries( gather PUBLIC core domain_decomposition field )
target_sources( gather
      PRIVATE
      CellGatherPackInfo.h
      FileGatherScheme.cpp
      GnuPlotGraphWriter.h
      CurveGatherPackInfo.h
      CurveGatherPackInfo.impl.h
      FileGatherScheme.h
      MPIGatherScheme.h
      CommonSchemeFunctions.h
      DataProcessor.h
      GatherPackInfo.h
      MPIGatherScheme.cpp
      CellGatherPackInfo.impl.h
      )
###################################################################################################
