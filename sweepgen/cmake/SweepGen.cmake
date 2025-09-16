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

#   CMake module for the waLBerla SweepGen metaprogramming system

set( 
    WALBERLA_SWEEPGEN_MANAGED_VENV ON
    CACHE BOOL 
    "Have SweepGen manage its virtual Python environment by itself" 
)
mark_as_advanced( WALBERLA_SWEEPGEN_MANAGED_VENV )

set( 
    WALBERLA_CODEGEN_SILENT_INSTALL ON
    CACHE BOOL 
    "Silence output of `pip install` during virtual environment setup" 
)
mark_as_advanced( WALBERLA_CODEGEN_SILENT_INSTALL )

if( WALBERLA_SWEEPGEN_MANAGED_VENV )
    find_package( Python COMPONENTS Interpreter REQUIRED )

    set(WALBERLA_CODEGEN_VENV_PATH ${walberla_BINARY_DIR}/sweepgen-venv CACHE PATH "Location of the virtual environment used for code generation")
    set(_venv_python_exe ${WALBERLA_CODEGEN_VENV_PATH}/bin/python)
    set(
        _WALBERLA_CODEGEN_VENV_STATEFILE
        ${CMAKE_CURRENT_BINARY_DIR}/walberla-venv-state.json
        CACHE INTERNAL
        "venv statefile - for internal use only"
    )
    list(APPEND
        _venv_manager_args
        -s ${_WALBERLA_CODEGEN_VENV_STATEFILE}
    )
    
    if( WALBERLA_CODEGEN_SILENT_INSTALL )
        list(APPEND _venv_manager_args -q)
    endif()

    set(
        _WALBERLA_CODEGEN_VENV_INVOKE_MANAGER
            ${Python_EXECUTABLE}
            ${walberla_SOURCE_DIR}/sweepgen/cmake/ManageCodegenVenv.py
            ${_venv_manager_args}
        CACHE INTERNAL
        "venv manager filepath - for internal use only"
    )

    if(DEFINED CACHE{SWEEPGEN_REQUIREMENTS_FILE})
        message(STATUS "Using custom requirements file for SweepGen: $CACHE{SWEEPGEN_REQUIREMENTS_FILE}")
    else()
        set(
            SWEEPGEN_REQUIREMENTS_FILE
            ${walberla_SOURCE_DIR}/sweepgen/cmake/sweepgen-requirements.txt
            CACHE PATH
            "Location of the primary requirements file for the SweepGen virtual environment"
        )
        mark_as_advanced(SWEEPGEN_REQUIREMENTS_FILE)
    endif()

    set(
        _requirements_file
        ${CMAKE_CURRENT_BINARY_DIR}/sweepgen-requirements.txt
    )

    file(COPY_FILE ${SWEEPGEN_REQUIREMENTS_FILE} ${_requirements_file} )
    file(APPEND ${_requirements_file} "\n-e ${walberla_SOURCE_DIR}/sweepgen\n")

    set(_init_args ${WALBERLA_CODEGEN_VENV_PATH} ${_requirements_file})

    if(NOT _walberla_codegen_venv_initialized)
        #   Force-rebuild environment if cmake cache was deleted
        list(APPEND _init_args --force)
    endif()

    message(STATUS "Preparing SweepGen Python environment (this might take a while)")

    execute_process(
        COMMAND
            $CACHE{_WALBERLA_CODEGEN_VENV_INVOKE_MANAGER}
            init
            ${_init_args}
        RESULT_VARIABLE _lastResult
        ERROR_VARIABLE _lastError
    )

    if( ${_lastResult} )
        message( FATAL_ERROR "Codegen virtual environment setup failed:\n${_lastError}" )
    endif()

    set(
        _walberla_codegen_venv_initialized
        TRUE
        CACHE INTERNAL "venv initialized - for internal use only"
    )

    set( _wlb_codegen_python_init ${_venv_python_exe} )
else()
    #   Use the external Python environment, but check if all packages are installed
    find_package( Python COMPONENTS Interpreter REQUIRED )

    execute_process(
        COMMAND ${Python_EXECUTABLE} ${walberla_SOURCE_DIR}/sweepgen/cmake/CheckExternalEnv.py
        RESULT_VARIABLE _package_check_status
        OUTPUT_VARIABLE _package_check_output
        ERROR_VARIABLE _package_check_error
    )

    message( STATUS ${_package_check_output} )

    if(NOT (${_package_check_status} EQUAL 0) )
        message( FATAL_ERROR ${_package_check_error} )
    endif()

    set( _wlb_codegen_python_init ${Python_EXECUTABLE} )
endif()

set(PystencilsSfg_PYTHON_PATH ${_wlb_codegen_python_init})
set(
    WALBERLA_CODEGEN_PYTHON
    ${_wlb_codegen_python_init}
    CACHE PATH
    "Path to the Python interpreter used for code generation"
)
mark_as_advanced(WALBERLA_CODEGEN_PYTHON)

#   Find pystencils-sfg

find_package( PystencilsSfg REQUIRED )

#   Project Configuration Module

set( 
    WALBERLA_CODEGEN_CONFIG_MODULE
    ${CMAKE_BINARY_DIR}/SweepgenConfig.py
    CACHE
    FILEPATH
    "Path to waLBerla-wide codegen config module" 
)
mark_as_advanced( WALBERLA_CODEGEN_CONFIG_MODULE )

configure_file(
    ${walberla_SOURCE_DIR}/sweepgen/cmake/SweepgenConfig.template.py
    ${WALBERLA_CODEGEN_CONFIG_MODULE}
)

message( STATUS "Wrote project-wide code generator configuration to ${WALBERLA_CODEGEN_CONFIG_MODULE}" )

function(sweepgen_venv_require)
    if(NOT WALBERLA_SWEEPGEN_MANAGED_VENV)
        return()
    endif()

    execute_process(
        COMMAND
            $CACHE{_WALBERLA_CODEGEN_VENV_INVOKE_MANAGER}
            require
            --
            ${ARGV}
        RESULT_VARIABLE _lastResult
        ERROR_VARIABLE _lastError
    )

    if( ${_lastResult} )
        message( FATAL_ERROR ${_lastError} )
    endif()
endfunction()

function(sweepgen_venv_populate)
    if(NOT WALBERLA_SWEEPGEN_MANAGED_VENV)
        return()
    endif()

    execute_process(
        COMMAND
            $CACHE{_WALBERLA_CODEGEN_VENV_INVOKE_MANAGER}
            populate
        RESULT_VARIABLE _lastResult
        ERROR_VARIABLE _lastError
    )

    if( ${_lastResult} )
        message( FATAL_ERROR ${_lastError} )
    endif()
endfunction()

#   Code Generation Functions

#[[
Register code generation scripts for a CMake target.

Signature:

```
walberla_generate_sources( <target> 
    SCRIPTS script1.py [script2.py ...]
    [DEPENDS dependency1.py [dependency2.py...] ]
    [FILE_EXTENSIONS <header-extension> <impl-extension>]
    [OUTPUT_MODE <standalone|inline|header-only>]
)
```

This is a wrapper around `pystencilssfg_generate_target_sources`
without the `CONFIG_MODULE` parameter.
See also https://pycodegen.pages.i10git.cs.fau.de/pystencils-sfg/usage/project_integration.html#add-generator-scripts
#]]
function(walberla_generate_sources TARGET)
    pystencilssfg_generate_target_sources(${ARGV} CONFIG_MODULE $CACHE{WALBERLA_CODEGEN_CONFIG_MODULE})
    target_link_libraries( ${TARGET} PUBLIC walberla_sweepgen )
endfunction()
