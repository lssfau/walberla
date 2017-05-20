waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

waLBerla_add_executable ( NAME 01_GameOfLife_cuda
                          FILES 01_GameOfLife_cuda.cpp 01_GameOfLife_kernels.cu
                          DEPENDS blockforest core cuda field lbm geometry timeloop gui )                          
                                