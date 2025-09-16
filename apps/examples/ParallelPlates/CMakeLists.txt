walberla_link_files_to_builddir( params )

set( AppName ParallelPlates )

add_executable( ${AppName} )
target_sources( ${AppName} PRIVATE ${AppName}.cpp )

walberla_generate_sources(
    ${AppName}
    SCRIPTS ParallelPlatesSweeps.py
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

add_dependencies( waLBerlaTestsuite ${AppName} )

walberla_execute_test(
    NO_MODULE_LABEL
    NAME "ParallelPlates - CPU - Couette"
    COMMAND $<TARGET_FILE:${AppName}> params/Couette.prm 
    PROCESSES 1
)

walberla_execute_test(
    NO_MODULE_LABEL
    NAME "ParallelPlates - CPU - Poiseuille"
    COMMAND $<TARGET_FILE:${AppName}> params/Poiseuille.prm 
    PROCESSES 1
)
