waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable ( NAME PE_GranularGas
                          FILES PE_GranularGas.cpp
                          DEPENDS blockforest core pe postprocessing )

waLBerla_add_executable ( NAME MESA_PD_GranularGas
                          FILES MESA_PD_GranularGas.cpp
                          DEPENDS blockforest core pe mesa_pd postprocessing vtk )

waLBerla_add_executable ( NAME MESA_PD_KernelBenchmark
                          FILES MESA_PD_KernelBenchmark.cpp
                          DEPENDS blockforest core pe mesa_pd postprocessing vtk )
