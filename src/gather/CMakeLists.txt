###################################################################################################
#
# Module gather ( all-to-one walberla::communication    for output )
#
###################################################################################################

add_library( walberla_gather )
add_library( walberla::gather ALIAS walberla_gather )
target_link_libraries( walberla_gather PUBLIC walberla::core walberla::domain_decomposition walberla::field )
target_sources( walberla_gather
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
