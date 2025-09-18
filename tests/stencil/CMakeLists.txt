###################################################################################################
#
# Tests for walberla::stencil module
#
###################################################################################################


waLBerla_add_test_executable ( StencilTest StencilTest.cpp )
target_link_libraries ( StencilTest walberla::stencil )
waLBerla_execute_test( NAME StencilTest )
