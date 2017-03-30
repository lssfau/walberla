
###################################################################################################
#
# Module physics engine (pe)
#
###################################################################################################

if(WALBERLA_DOUBLE_ACCURACY)
    set(CCD_DOUBLE ON)
    set(CCD_SINGLE OFF)
else()
    set(CCD_DOUBLE OFF)
    set(CCD_SINGLE ON)
endif()

configure_file(extern/libccd/ccd/config.h.cmake.in extern/libccd/ccd/config.h)

include_directories(extern/libccd)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/extern/libccd")

waLBerla_add_module( DEPENDS core blockforest domain_decomposition geometry stencil vtk  )

if( WALBERLA_CXX_COMPILER_IS_MSVC )
   SET_SOURCE_FILES_PROPERTIES( extern/libccd/ccd.c
                                extern/libccd/mpr.c
                                extern/libccd/polytope.c
                                extern/libccd/support.c
                                extern/libccd/vec3.c
                                PROPERTIES LANGUAGE CXX )
endif()

