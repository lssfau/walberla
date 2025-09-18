###################################################################################################
#
# Tests for geometry_input module
#
###################################################################################################

# Test for SendBuffer and RecvBuffer
waLBerla_add_test_executable ( CellAABBTest CellAABBTest.cpp )
target_link_libraries ( CellAABBTest walberla::geometry )
waLBerla_execute_test( NAME CellAABBTest COMMAND $<TARGET_FILE:CellAABBTest> )

waLBerla_add_test_executable ( Functions Functions.cpp )
target_link_libraries ( Functions walberla::geometry )
waLBerla_execute_test( NAME Functions )

waLBerla_add_test_executable ( VoxelFileTest VoxelFileTest.cpp )
target_link_libraries ( VoxelFileTest walberla::geometry )
waLBerla_execute_test( NAME VoxelFileTest     COMMAND $<TARGET_FILE:VoxelFileTest> )
waLBerla_execute_test( NAME VoxelFileTestLong COMMAND $<TARGET_FILE:VoxelFileTest> --longrun LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )
set_property( TEST VoxelFileTestLong PROPERTY DEPENDS VoxelFileTest ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN

waLBerla_add_test_executable ( ScalarFieldFromBodyTest ScalarFieldFromBodyTest.cpp )
target_link_libraries ( ScalarFieldFromBodyTest walberla::geometry )
waLBerla_execute_test( NAME ScalarFieldFromBodyTest )


file( COPY "test.png" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
waLBerla_add_test_executable ( ScalarFieldFromGrayScaleImageTest ScalarFieldFromGrayScaleImageTest.cpp )
target_link_libraries ( ScalarFieldFromGrayScaleImageTest walberla::geometry )
waLBerla_execute_test( NAME ScalarFieldFromGrayScaleImageTest )


waLBerla_add_test_executable ( ImageTest ImageTest.cpp )
target_link_libraries ( ImageTest walberla::geometry )
waLBerla_execute_test( NAME ImageTest )

waLBerla_add_test_executable ( BinaryRawFileTest BinaryRawFileTest.cpp )
target_link_libraries ( BinaryRawFileTest walberla::geometry )
