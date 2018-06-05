waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable ( NAME PeriodicGranularGas
                          FILES PeriodicGranularGas.cpp
                          DEPENDS blockforest core pe )
