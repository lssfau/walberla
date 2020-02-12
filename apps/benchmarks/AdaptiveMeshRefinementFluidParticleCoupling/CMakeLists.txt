waLBerla_add_executable( NAME WorkloadEvaluation FILES WorkloadEvaluation.cpp DEPENDS blockforest boundary core field lbm pe pe_coupling postprocessing stencil timeloop vtk )

waLBerla_add_executable( NAME AMRSedimentSettling FILES AMRSedimentSettling.cpp DEPENDS blockforest boundary core field lbm pe pe_coupling postprocessing stencil timeloop vtk )

waLBerla_add_executable( NAME AMRSettlingSphere FILES AMRSettlingSphere.cpp DEPENDS blockforest boundary core field lbm pe pe_coupling postprocessing stencil timeloop vtk )
