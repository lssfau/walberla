waLBerla_add_executable( NAME WorkloadEvaluation FILES WorkloadEvaluation.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm walberla::pe walberla::pe_coupling walberla::postprocessing walberla::stencil walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME AMRSedimentSettling FILES AMRSedimentSettling.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm walberla::pe walberla::pe_coupling walberla::postprocessing walberla::stencil walberla::timeloop walberla::vtk )

waLBerla_add_executable( NAME AMRSettlingSphere FILES AMRSettlingSphere.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm walberla::pe walberla::pe_coupling walberla::postprocessing walberla::stencil walberla::timeloop walberla::vtk )
