# waLBerla Python module


if ( WALBERLA_BUILD_WITH_PYTHON_MODULE )
    
    if( WALBERLA_CXX_COMPILER_IS_MSVC )
       set ( pythonModules blockforest boundary domain_decomposition core field geometry lbm postprocessing python_coupling timeloop vtk)
    elseif( APPLE )
       set ( pythonModules "-Wl,-force_load" blockforest boundary domain_decomposition core field geometry lbm postprocessing python_coupling timeloop vtk)
    else()
       set ( pythonModules "-Wl,-whole-archive" blockforest boundary domain_decomposition core field geometry lbm postprocessing python_coupling timeloop vtk  "-Wl,-no-whole-archive" )
    endif()


    #add_library( walberla_cpp SHARED PythonModuleWithLbmModule.cpp )
    add_library( walberla_cpp SHARED PythonModule.cpp )

    target_link_libraries( walberla_cpp ${pythonModules} ${SERVICE_LIBS}  )
    set_target_properties( walberla_cpp PROPERTIES PREFIX "")
    if ( APPLE )
       set_target_properties( walberla_cpp PROPERTIES SUFFIX ".so")
    endif()
    set_target_properties( walberla_cpp PROPERTIES MACOSX_RPATH TRUE )

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/setup.py  ${CMAKE_CURRENT_BINARY_DIR}/setup.py )

    set( PYTHON_EXE /usr/bin/env python)
    if( ${PYTHON_LIBRARY} MATCHES "python3" )
        set( PYTHON_EXE /usr/bin/env python3 )
    endif()

    add_custom_target( pythonModule        ALL ${PYTHON_EXE} setup.py build  DEPENDS walberla_cpp  )
    add_custom_target( pythonModuleInstall     ${PYTHON_EXE} setup.py install DEPENDS walberla_cpp )

    add_test( NAME PythonModuleTest
              COMMAND ${PYTHON_EXE} -m unittest discover -v -s ${walberla_SOURCE_DIR}/python/waLBerla_tests )
endif()