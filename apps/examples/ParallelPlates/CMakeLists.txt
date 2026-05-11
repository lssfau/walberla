walberla_link_files_to_builddir( params )

set( AppName ParallelPlates )

add_executable( ${AppName} )
target_sources( ${AppName} PRIVATE ${AppName}.cpp )
walberla_set_gpu_language( ${AppName}.cpp )

walberla_generate_sources(
    ${AppName}
    SCRIPTS ParallelPlatesSweeps.py
    AUTO_LANGUAGE
)

target_link_libraries(
    ${AppName}
    PRIVATE
    walberla::v8
)

add_dependencies( waLBerlaTestsuite ${AppName} )

walberla_execute_test(
    NO_MODULE_LABEL
    NAME "ParallelPlates - Couette"
    COMMAND $<TARGET_FILE:${AppName}> params/Couette.prm 
    PROCESSES 1
)

walberla_execute_test(
    NO_MODULE_LABEL
    NAME "ParallelPlates - Poiseuille"
    COMMAND $<TARGET_FILE:${AppName}> params/Poiseuille.prm 
    PROCESSES 1
)
