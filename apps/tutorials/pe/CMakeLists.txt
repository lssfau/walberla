waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable ( NAME 01_Tutorial_ConfinedGas 
                          FILES 01_ConfinedGas.cpp
                          DEPENDS blockforest core pe )
