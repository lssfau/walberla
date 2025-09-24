waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

add_executable ( 01_GameOfLife_cuda 01_GameOfLife_cuda.cpp 01_GameOfLife_kernels.cu )
target_link_libraries( 01_GameOfLife_cuda walberla::blockforest walberla::core walberla::gpu walberla::field
        walberla::lbm walberla::geometry walberla::timeloop )
