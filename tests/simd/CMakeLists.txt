###################################################################################################
#
# Tests for simd module
#
###################################################################################################



# builds the test

if ( CMAKE_COMPILER_IS_GNUCXX )
   if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0" )
      set( MarchNativeString "" )
   else()
      set( MarchNativeString "-march=native" )
   endif()
endif()

waLBerla_compile_test( NAME   AVX2_AVX_Equivalence FILES SIMD_Equivalence.cpp  )
set_property         ( TARGET AVX2_AVX_Equivalence PROPERTY COMPILE_FLAGS "-DIS0_AVX2 -DIS1_AVX" )
waLBerla_execute_test( NAME   AVX2_AVX_Equivalence )


waLBerla_compile_test( NAME   AVX_SSE4_Equivalence  FILES SIMD_Equivalence.cpp  )
set_property         ( TARGET AVX_SSE4_Equivalence  PROPERTY COMPILE_FLAGS "${MarchNativeString} -DIS0_AVX -DIS1_SSE4" )
waLBerla_execute_test( NAME   AVX_SSE4_Equivalence )


waLBerla_compile_test( NAME   SSE4_SSE2_Equivalence FILES SIMD_Equivalence.cpp  )
set_property         ( TARGET SSE4_SSE2_Equivalence PROPERTY COMPILE_FLAGS "${MarchNativeString} -DIS0_SSE2 -DIS1_SSE4" )
waLBerla_execute_test( NAME   SSE4_SSE2_Equivalence )


waLBerla_compile_test( NAME   AVX_Scalar_Equivalence FILES SIMD_Equivalence.cpp  )
set_property         ( TARGET AVX_Scalar_Equivalence PROPERTY COMPILE_FLAGS "${MarchNativeString} -DIS0_AVX -DIS1_SCALAR" )
waLBerla_execute_test( NAME   AVX_Scalar_Equivalence )


waLBerla_compile_test( NAME   QPX_Scalar_Equivalence FILE SIMD_Equivalence.cpp  )
set_property         ( TARGET QPX_Scalar_Equivalence PROPERTY COMPILE_FLAGS "-DIS0_QPX -DIS1_SCALAR" )
waLBerla_execute_test( NAME   QPX_Scalar_Equivalence )

