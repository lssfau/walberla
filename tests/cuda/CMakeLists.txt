###################################################################################################
#
# Tests for cuda
#
###################################################################################################

waLBerla_compile_test( FILES communication/GPUPackInfoTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME  GPUPackInfoTest )

waLBerla_compile_test( FILES communication/CommTest )
waLBerla_execute_test( NAME  CommTest )

waLBerla_compile_test( FILES FieldTransferTest )
waLBerla_execute_test( NAME  FieldTransferTest )

waLBerla_compile_test( FILES SimpleKernelTest.cpp Kernels.cu DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME  SimpleKernelTest )

waLBerla_compile_test( FILES CudaMPI DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME  CudaMPI )

waLBerla_compile_test( FILES FieldIndexing3DTest.cpp FieldIndexing3DTest.cu )
waLBerla_execute_test( NAME  FieldIndexing3DTest )

