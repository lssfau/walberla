waLBerla_link_files_to_builddir( *.prm )

add_executable ( DeformationField DeformationField.cpp )
target_link_libraries ( DeformationField walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )

add_executable ( SingleVortex SingleVortex.cpp )
target_link_libraries ( SingleVortex walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing
        walberla::timeloop walberla::vtk )

add_executable ( ZalesakDisk ZalesakDisk.cpp )
target_link_libraries ( ZalesakDisk walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk )