# This file is part of waLBerla. waLBerla is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# waLBerla is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
#
# Author: Frederik Hennig <frederik.hennig@fau.de>
#

#[[

Targets and Build Settings for the waLBerla v8 Core Library
===========================================================

]]#

add_library( walberla_v8 INTERFACE )
add_library( walberla::v8 ALIAS walberla_v8 )

target_include_directories(
    walberla_v8
    INTERFACE
    ${walberla_SOURCE_DIR}/include
)

target_link_libraries(
    walberla_v8
    INTERFACE
    walberla_core
    walberla_stencil
    walberla_blockforest
    walberla_field
)

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_HIP )
    target_link_libraries( walberla_v8 INTERFACE walberla_gpu )
endif()

if( WALBERLA_BUILD_WITH_CUDA )

    #   Enable relaxed constexpr for consuming CUDA TUs
    target_compile_options(
        walberla_v8
        INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    )

elseif( WALBERLA_BUILD_WITH_HIP )
    
    #   Fetch and link core against HIP C++ standard libraries
    FetchContent_Declare(
        libhipcxx
        GIT_REPOSITORY https://github.com/ROCm/libhipcxx.git
        GIT_TAG release/2.7.x
    )

    FetchContent_MakeAvailable( libhipcxx ) 

    target_link_libraries( walberla_core PUBLIC libhipcxx::libhipcxx )

    #   Suppress custom `#warning` messages in hip/std runtime headers which would otherwise cause compilation to fail
    #   See https://clang.llvm.org/docs/WarningSuppressionMappings.html
    target_compile_options(
        walberla_v8
        INTERFACE
        $<$<COMPILE_LANGUAGE:HIP>:--warning-suppression-mappings=${walberla_SOURCE_DIR}/cmake/hip-warning-suppression.txt>
    )

endif()

#######################################################################################################################
#
# Set the compilation language of source files to the current waLBerla GPU language,
# and modifies the compilation options of the given source files to accomodate the portable waLBerla runtime library.
#
#######################################################################################################################

function(walberla_set_gpu_language FILES)
    if( $CACHE{WALBERLA_BUILD_WITH_CUDA} )
        set_source_files_properties(
            ${FILES}
            PROPERTIES 
            LANGUAGE CUDA
        )
    elseif( $CACHE{WALBERLA_BUILD_WITH_HIP} )
        set_source_files_properties(
            ${FILES}
            PROPERTIES
            LANGUAGE HIP
        )
    endif()
endfunction()
