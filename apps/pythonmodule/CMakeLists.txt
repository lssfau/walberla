# waLBerla Python module
if( NOT TARGET python_coupling )
    message( WARNING "python module ist not build since the python_coupling target is non-existent" )
else()
if ( WALBERLA_BUILD_WITH_PYTHON )
    if ( WALBERLA_BUILD_WITH_GPU_SUPPORT )
        set(PYTHON_MODULE_DEPENDENCIES blockforest boundary domain_decomposition core field python_coupling timeloop vtk gpu)
    else()
        set(PYTHON_MODULE_DEPENDENCIES blockforest boundary domain_decomposition core field python_coupling timeloop vtk)
    endif()

    if( WALBERLA_CXX_COMPILER_IS_MSVC )
        set ( pythonModules ${PYTHON_MODULE_DEPENDENCIES})
    elseif( APPLE )
        set ( pythonModules "-Wl,-force_load" ${PYTHON_MODULE_DEPENDENCIES})
    else()
        set ( pythonModules "-Wl,-whole-archive" ${PYTHON_MODULE_DEPENDENCIES}  "-Wl,-no-whole-archive" )
    endif()

    add_library( walberla_cpp SHARED PythonModule.cpp )

    target_link_libraries( walberla_cpp ${WALBERLA_LINK_LIBRARIES_KEYWORD} ${pythonModules} ${SERVICE_LIBS} )

    set_target_properties(
            walberla_cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${PYTHON_MODULE_EXTENSION}"
    )
    set_target_properties( walberla_cpp PROPERTIES MACOSX_RPATH TRUE )

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/setup.py  ${CMAKE_CURRENT_BINARY_DIR}/setup.py )

    add_custom_target( pythonModule        ALL ${Python_EXECUTABLE} setup.py build   DEPENDS walberla_cpp )
    add_custom_target( pythonModuleInstall     ${Python_EXECUTABLE} setup.py install --user DEPENDS walberla_cpp )
endif()
endif()