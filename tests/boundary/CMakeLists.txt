###################################################################################################
#
# Tests for boundary module
#
###################################################################################################


waLBerla_compile_test( FILES BoundaryHandling.cpp )
waLBerla_execute_test( NAME BoundaryHandling )
