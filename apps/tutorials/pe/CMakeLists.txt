waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable ( NAME 01_Tutorial_ConfinedGas
                          FILES 01_ConfinedGas.cpp
                          DEPENDS blockforest core pe )

waLBerla_add_executable ( NAME 02_Tutorial_ConfinedGasExtended
                          FILES 02_ConfinedGasExtended.cpp
                          DEPENDS blockforest core pe postprocessing vtk sqlite )

waLBerla_execute_test( NO_MODULE_LABEL NAME 01_Tutorial_ConfinedGas PROCESSES 8 )
waLBerla_execute_test( NO_MODULE_LABEL NAME 02_Tutorial_ConfinedGasExtended
                       COMMAND $<TARGET_FILE:02_Tutorial_ConfinedGasExtended> 02_ConfinedGasExtended.cfg
                       PROCESSES 8 DEPENDS_ON_TARGETS 02_Tutorial_ConfinedGasExtended )
