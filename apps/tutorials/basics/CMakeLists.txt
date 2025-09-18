waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

add_executable ( 01_Tutorial_BlocksAndFields1 01_BlocksAndFields.cpp )
target_link_libraries ( 01_Tutorial_BlocksAndFields1 walberla::blockforest walberla::core walberla::field
        walberla::timeloop )
                                                                              
add_executable ( 02_Tutorial_Sweeps 02_Sweeps.cpp )
target_link_libraries ( 02_Tutorial_Sweeps walberla::blockforest walberla::core walberla::field walberla::timeloop )
  
add_executable ( 03_Tutorial_GameOfLife 03_GameOfLife.cpp )
target_link_libraries ( 03_Tutorial_GameOfLife walberla::blockforest walberla::core walberla::field walberla::geometry
        walberla::timeloop )
                                                                                             
