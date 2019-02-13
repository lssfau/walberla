# waLBerla Python module


if ( WALBERLA_BUILD_WITH_PYTHON_MODULE )

    set(PYTHON_MODULE_DEPENDENCIES blockforest boundary domain_decomposition core field geometry lbm postprocessing python_coupling timeloop vtk)
    if (WALBERLA_BUILD_WITH_CUDA)
        set(PYTHON_MODULE_DEPENDENCIES ${PYTHON_MODULE_DEPENDENCIES} cuda)
    endif()

    if (WALBERLA_BUILD_WITH_OPENMESH)
        set(PYTHON_MODULE_DEPENDENCIES ${PYTHON_MODULE_DEPENDENCIES} mesh pe)
    endif()

    if( WALBERLA_CXX_COMPILER_IS_MSVC )
       set ( pythonModules ${PYTHON_MODULE_DEPENDENCIES})
    elseif( APPLE )
       set ( pythonModules "-Wl,-force_load" ${PYTHON_MODULE_DEPENDENCIES})
    else()
       set ( pythonModules "-Wl,-whole-archive" ${PYTHON_MODULE_DEPENDENCIES}  "-Wl,-no-whole-archive" )
    endif()

    if( WALBERLA_BUILD_WITH_PYTHON_LBM )
        add_library( walberla_cpp SHARED PythonModuleWithLbmModule.cpp )
    else()
        add_library( walberla_cpp SHARED PythonModule.cpp )
    endif()

    target_link_libraries( walberla_cpp ${WALBERLA_LINK_LIBRARIES_KEYWORD} ${pythonModules} ${SERVICE_LIBS}  )
    set_target_properties( walberla_cpp PROPERTIES PREFIX "")
    if ( APPLE )
       set_target_properties( walberla_cpp PROPERTIES SUFFIX ".so")
    endif()
    set_target_properties( walberla_cpp PROPERTIES MACOSX_RPATH TRUE )

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/setup.py  ${CMAKE_CURRENT_BINARY_DIR}/setup.py )

    add_custom_target( pythonModule        ALL ${PYTHON_EXECUTABLE} setup.py build   DEPENDS walberla_cpp )
    add_custom_target( pythonModuleInstall     ${PYTHON_EXECUTABLE} setup.py install DEPENDS walberla_cpp )

    add_test( NAME PythonModuleTest
              COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -v -s ${walberla_SOURCE_DIR}/python/waLBerla_tests )
endif()
