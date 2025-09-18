add_library( walberla_fft )
add_library( walberla::fft ALIAS walberla_fft )
target_link_libraries( walberla_fft PUBLIC walberla::core walberla::blockforest walberla::field )

if ( PFFT_FOUND AND FFTW3_MPI_FOUND )
    target_link_libraries( walberla_fft PUBLIC ${PFFT_LIBRARIES} ${FFTW3_LIBRARIES} ${FFTW3_MPI_LIBRARIES} )
elseif ( FFTW3_FOUND )
    target_link_libraries( walberla_fft PUBLIC ${FFTW3_LIBRARIES} )
endif ()

target_sources( walberla_fft PUBLIC Fft.h )
target_sources( walberla_fft PRIVATE Fft.cpp )

