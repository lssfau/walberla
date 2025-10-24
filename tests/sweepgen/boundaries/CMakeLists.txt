walberla_link_files_to_builddir( *.prm )
waLBerla_link_geometry_to_builddir(cube.obj)

set( AppName TestBoundaryLinks )

waLBerla_add_test_executable( ${AppName} )
target_sources( ${AppName} PRIVATE ${AppName}.cpp )

walberla_generate_sources(
    ${AppName}
    SCRIPTS ${AppName}.py
)

target_link_libraries(
    ${AppName}
    PRIVATE
    walberla::core
    walberla::experimental
    walberla::blockforest
    walberla::geometry
    walberla::stencil
)

if( WALBERLA_BUILD_WITH_OPENMESH )
    target_link_libraries(
            ${AppName}
            PRIVATE
            walberla::mesh
            walberla::mesh_common
    )
endif()

waLBerla_execute_test(NAME ${AppName} COMMAND ${AppName})