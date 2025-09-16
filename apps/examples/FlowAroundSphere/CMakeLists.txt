walberla_link_files_to_builddir( *.prm )

set( AppName FlowAroundSphereExample )

add_executable( ${AppName} )
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
    walberla::vtk
)
