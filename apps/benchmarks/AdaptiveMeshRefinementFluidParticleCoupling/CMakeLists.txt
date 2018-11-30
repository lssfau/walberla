waLBerla_add_executable( NAME WorkloadEvaluation FILES WorkloadEvaluation.cpp DEPENDS blockforest boundary core field lbm pe pe_coupling postprocessing stencil timeloop vtk )

if( WALBERLA_BUILD_WITH_PARMETIS )
   waLBerla_add_executable( NAME AMRSedimentSettling FILES AMRSedimentSettling.cpp DEPENDS blockforest boundary core field lbm pe pe_coupling postprocessing stencil timeloop vtk )
endif()