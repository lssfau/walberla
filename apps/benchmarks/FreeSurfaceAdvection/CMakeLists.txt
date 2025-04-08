waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable(NAME    DeformationField
                        FILES   DeformationField.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable(NAME    SingleVortex
                        FILES   SingleVortex.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )

waLBerla_add_executable(NAME    ZalesakDisk
                        FILES   ZalesakDisk.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk )