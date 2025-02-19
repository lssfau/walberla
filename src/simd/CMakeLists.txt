
###################################################################################################
#
# Module phasefield
#
###################################################################################################

add_library( walberla_simd INTERFACE )
add_library( walberla::simd ALIAS walberla_simd )
target_link_libraries( walberla_simd INTERFACE walberla::core )
target_sources( walberla_simd
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

