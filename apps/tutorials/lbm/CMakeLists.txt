waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.png )

add_executable ( 01_BasicLBM 01_BasicLBM.cpp )
target_link_libraries ( 01_BasicLBM walberla::blockforest walberla::core walberla::field walberla::lbm
        walberla::geometry walberla::timeloop )

add_executable ( 02_BasicLBM_ExemplaryExtensions 02_BasicLBM_ExemplaryExtensions.cpp )
target_link_libraries ( 02_BasicLBM_ExemplaryExtensions walberla::blockforest walberla::core walberla::field
        walberla::lbm walberla::geometry walberla::timeloop )

add_executable ( 03_LBLidDrivenCavity 03_LBLidDrivenCavity.cpp )
target_link_libraries ( 03_LBLidDrivenCavity walberla::blockforest walberla::boundary walberla::core walberla::field
        walberla::lbm walberla::stencil walberla::timeloop walberla::vtk )

if ( WALBERLA_BUILD_WITH_OPENMESH )

waLBerla_link_geometry_to_builddir(*.obj)

add_executable ( 04_LBComplexGeometry 04_LBComplexGeometry.cpp )
target_link_libraries( 04_LBComplexGeometry walberla::blockforest walberla::boundary walberla::core walberla::field
        walberla::lbm walberla::mesh walberla::stencil walberla::timeloop walberla::vtk )

endif()

add_executable ( 05_BackwardFacingStep 05_BackwardFacingStep.cpp )
target_link_libraries ( 05_BackwardFacingStep walberla::blockforest walberla::boundary walberla::core walberla::field
        walberla::lbm walberla::stencil walberla::timeloop walberla::vtk )

add_executable ( 06_LBBoundaryCondition 06_LBBoundaryCondition.cpp )
target_link_libraries ( 06_LBBoundaryCondition walberla::blockforest walberla::boundary walberla::core walberla::field
        walberla::lbm walberla::stencil walberla::timeloop walberla::vtk )