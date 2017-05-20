###################################################################################################
#
# Tests for cuda
#
###################################################################################################

waLBerla_compile_test( FILES FieldTransferTest )
waLBerla_execute_test( NAME  FieldTransferTest )

waLBerla_compile_test( FILES SimpleKernelTest.cpp Kernels.cu DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME  SimpleKernelTest )

waLBerla_compile_test( FILES CudaMPI DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME  CudaMPI )

