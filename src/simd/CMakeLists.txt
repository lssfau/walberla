
###################################################################################################
#
# Module phasefield
#
###################################################################################################

add_library( simd INTERFACE )
target_link_libraries( simd INTERFACE core )
target_sources( simd
      INTERFACE
      AVX.h
      SSE4.h
      Scalar.h
      AlignedAllocator.h
      StreamOutput.h
      IntelVecTypesCppOperators.h
      SIMD.h
      AVX2.h
      SSE2.h
      )

