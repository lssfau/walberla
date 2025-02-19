waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

waLBerla_add_executable ( NAME 01_GameOfLife_cuda
                          FILES 01_GameOfLife_cuda.cpp 01_GameOfLife_kernels.cu
      DEPENDS walberla::blockforest walberla::core walberla::gpu walberla::field walberla::lbm walberla::geometry walberla::timeloop )
                                