waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )


waLBerla_add_executable ( NAME 01_SolvingPDE 
                          FILES 01_SolvingPDE.cpp
                          DEPENDS blockforest core field stencil timeloop gui )
                          
waLBerla_add_executable ( NAME 02_HeatEquation 
                          FILES 02_HeatEquation.cpp
                          DEPENDS blockforest core field stencil timeloop gui )
                          
waLBerla_add_executable ( NAME 03_HeatEquation_Extensions 
                          FILES 03_HeatEquation_Extensions.cpp
                          DEPENDS blockforest core field stencil timeloop gui vtk )
                          