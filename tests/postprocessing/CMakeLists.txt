###################################################################################################
#
# Tests for walberla::postprocessing module
#
###################################################################################################


waLBerla_add_test_executable ( SphereTriangulate SphereTriangulate.cpp )
target_link_libraries ( SphereTriangulate walberla::postprocessing )
waLBerla_execute_test( NAME SphereTriangulate )