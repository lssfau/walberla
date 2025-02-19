###################################################################################################
#
# Module fft 
#
###################################################################################################

add_library( walberla_fft )
add_library( walberla::fft ALIAS walberla_fft )
target_link_libraries( walberla_fft PUBLIC walberla::core walberla::blockforest walberla::field )
target_sources( walberla_fft
      PRIVATE
      Fft.h
      Fft.cpp
      )

###################################################################################################
