###################################################################################################
#
# Tests for postprocessing module
#
###################################################################################################


waLBerla_compile_test( FILES SphereTriangulate.cpp )
waLBerla_execute_test( NAME SphereTriangulate )

waLBerla_compile_test( FILES SQLiteTest.cpp )
waLBerla_execute_test( NAME SQLiteTest )
