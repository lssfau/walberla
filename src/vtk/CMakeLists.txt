###################################################################################################
#
# Module vtk
#
###################################################################################################

add_library( vtk )
target_sources( vtk
      PRIVATE
      Base64Writer.cpp
      Base64Writer.h
      DumpBlockStructureLevel.h
      CellBBCellFilter.h
      Initialization.h
      Initialization.cpp
      PointDataSource.h
      BlockCellDataWriter.h
      VTKTrait.h
      VTKOutput.cpp
      VTKOutput.h
      DumpBlockStructureProcess.h
      AABBCellFilter.h
      PolylineDataSource.h
      UtilityFunctions.h
      ChainedFilter.h
      )
target_link_libraries( vtk PUBLIC core blockforest domain_decomposition )

###################################################################################################
