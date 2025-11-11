# waLBerla Python module
if ( NOT TARGET walberla::python_coupling )
    message( WARNING "python module ist not build since the walberla::python_coupling target is non-existent" )
else()
if ( WALBERLA_BUILD_WITH_PYTHON )
    if ( WALBERLA_BUILD_WITH_GPU_SUPPORT )
        set( PYTHON_MODULE_DEPENDENCIES walberla::blockforest walberla::boundary walberla::domain_decomposition walberla::core walberla::field walberla::python_coupling walberla::timeloop walberla::vtk walberla::gpu )
    else()
        set( PYTHON_MODULE_DEPENDENCIES walberla::blockforest walberla::boundary walberla::domain_decomposition walberla::core walberla::field walberla::python_coupling walberla::timeloop walberla::vtk )
    endif()

    if( APPLE )
        set ( pythonModules "-Wl,-force_load" ${PYTHON_MODULE_DEPENDENCIES})
    else()
        set ( pythonModules "-Wl,-whole-archive" ${PYTHON_MODULE_DEPENDENCIES}  "-Wl,-no-whole-archive" )
    endif()

    add_library( walberla_cpp SHARED PythonModule.cpp )

    target_link_libraries( walberla_cpp ${WALBERLA_LINK_LIBRARIES_KEYWORD} ${pythonModules} ${CORE_SERVICE_LIBS} )

    set_target_properties(
            walberla_cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${PYTHON_MODULE_EXTENSION}"
    )
    set_target_properties( walberla_cpp PROPERTIES MACOSX_RPATH TRUE )

    set(WALBERLA_MODULE_FILE_NAME "${PYTHON_MODULE_PREFIX}walberla_cpp${PYTHON_MODULE_EXTENSION}")

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml  ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml )
    add_custom_target( pythonModule ALL
                       COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
                               "${walberla_SOURCE_DIR}/python/waLBerla/"
                               "${CMAKE_CURRENT_BINARY_DIR}/waLBerla/"
                       COMMAND ${CMAKE_COMMAND} -E copy_if_different
                                "${CMAKE_CURRENT_BINARY_DIR}/${WALBERLA_MODULE_FILE_NAME}"
                                "${CMAKE_CURRENT_BINARY_DIR}/waLBerla/${WALBERLA_MODULE_FILE_NAME}"
                       DEPENDS walberla_cpp)

    add_custom_target( pythonModuleInstall ${Python_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR} DEPENDS pythonModule )
endif()
endif()