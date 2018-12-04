###################################################################################################
#
# Tests for geometry_input module
#
###################################################################################################

# Test for SendBuffer and RecvBuffer
waLBerla_compile_test( FILES CellAABBTest.cpp )
waLBerla_execute_test( NAME CellAABBTest COMMAND $<TARGET_FILE:CellAABBTest> )

waLBerla_compile_test( FILES Functions.cpp )
waLBerla_execute_test( NAME Functions )

waLBerla_compile_test( FILES VoxelFileTest.cpp )
waLBerla_execute_test( NAME VoxelFileTest     COMMAND $<TARGET_FILE:VoxelFileTest> )
waLBerla_execute_test( NAME VoxelFileTestLong COMMAND $<TARGET_FILE:VoxelFileTest> --longrun LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )
set_property( TEST VoxelFileTestLong PROPERTY DEPENDS VoxelFileTest ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN


waLBerla_compile_test( FILES ScalarFieldFromBodyTest.cpp DEPENDS gui )
waLBerla_execute_test( NAME ScalarFieldFromBodyTest )


file( COPY "test.png" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
waLBerla_compile_test( FILES ScalarFieldFromGrayScaleImageTest.cpp DEPENDS gui )
waLBerla_execute_test( NAME ScalarFieldFromGrayScaleImageTest )


waLBerla_compile_test( FILES ImageTest.cpp )
waLBerla_execute_test( NAME ImageTest )

waLBerla_compile_test( FILES BinaryRawFileTest.cpp )
