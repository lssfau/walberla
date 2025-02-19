###################################################################################################
#
# Module vtk
#
###################################################################################################

add_library( walberla_vtk )
add_library( walberla::vtk ALIAS walberla_vtk )
target_sources( walberla_vtk
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
target_link_libraries( walberla_vtk PUBLIC walberla::core walberla::blockforest walberla::domain_decomposition )

###################################################################################################
