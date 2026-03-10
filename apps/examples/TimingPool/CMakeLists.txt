walberla_link_files_to_builddir( *.prm )

set( AppName TimingPoolExample )

add_executable( ${AppName} )
target_sources( ${AppName} PRIVATE ${AppName}.cpp )

target_link_libraries(
    ${AppName}
    PRIVATE
    walberla::blockforest
    walberla::core
    walberla::field
    walberla::timeloop
)

add_dependencies( waLBerlaTestsuite ${AppName} )

walberla_execute_test(
    NO_MODULE_LABEL
    NAME "TimingPool - CPU - Serial"
    COMMAND $<TARGET_FILE:${AppName}> TimingPoolExample.prm
    PROCESSES 1
)
