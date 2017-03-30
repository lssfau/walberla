###################################################################################################
#
# Tests for stencil module
#
###################################################################################################


waLBerla_compile_test( FILES StencilTest.cpp )
waLBerla_execute_test( NAME StencilTest )
