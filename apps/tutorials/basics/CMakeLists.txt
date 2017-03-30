waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

waLBerla_add_executable ( NAME 01_Tutorial_BlocksAndFields1 
                          FILES 01_BlocksAndFields.cpp
                          DEPENDS blockforest core field timeloop gui )
                                                                              
waLBerla_add_executable ( NAME 02_Tutorial_Sweeps 
                          FILES 02_Sweeps.cpp
                          DEPENDS blockforest core field timeloop gui ) 
  
waLBerla_add_executable ( NAME 03_Tutorial_GameOfLife
                          FILES 03_GameOfLife.cpp
                          DEPENDS blockforest core field geometry timeloop gui )
                                                                                             
