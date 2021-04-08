waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

waLBerla_add_executable ( NAME 01_BasicLBM
                          FILES 01_BasicLBM.cpp
                          DEPENDS blockforest core field lbm geometry timeloop gui )
                                                                                             
waLBerla_add_executable ( NAME 02_BasicLBM_ExemplaryExtensions
                          FILES 02_BasicLBM_ExemplaryExtensions.cpp
                          DEPENDS blockforest core field lbm geometry timeloop gui )
                                                                                             
waLBerla_add_executable ( NAME 03_LBLidDrivenCavity 
                          FILES 03_LBLidDrivenCavity.cpp
                          DEPENDS blockforest boundary core field lbm stencil timeloop vtk )

if ( WALBERLA_BUILD_WITH_OPENMESH )

waLBerla_link_files_to_builddir( *.obj )

waLBerla_add_executable ( NAME 04_LBComplexGeometry
                          FILES 04_LBComplexGeometry.cpp
                          DEPENDS blockforest boundary core field lbm mesh stencil timeloop vtk )
                          
endif()

waLBerla_add_executable ( NAME 05_BackwardFacingStep
                          FILES 05_BackwardFacingStep.cpp
                          DEPENDS blockforest boundary core field lbm stencil timeloop vtk )

waLBerla_add_executable ( NAME 06_LBBoundaryCondition
                          FILES 06_LBBoundaryCondition.cpp
                          DEPENDS blockforest boundary core field lbm stencil timeloop vtk )