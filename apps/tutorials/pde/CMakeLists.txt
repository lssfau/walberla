waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )


add_executable ( 01_SolvingPDE 01_SolvingPDE.cpp )
target_link_libraries ( 01_SolvingPDE walberla::blockforest walberla::core walberla::field walberla::stencil
        walberla::timeloop )
                          
add_executable ( 02_HeatEquation 02_HeatEquation.cpp )
target_link_libraries( 02_HeatEquation walberla::blockforest walberla::core walberla::field walberla::stencil
        walberla::timeloop )
                          
add_executable ( 03_HeatEquation_Extensions 03_HeatEquation_Extensions.cpp )
target_link_libraries ( 03_HeatEquation_Extensions walberla::blockforest walberla::core walberla::field
        walberla::stencil walberla::timeloop walberla::vtk )
                          