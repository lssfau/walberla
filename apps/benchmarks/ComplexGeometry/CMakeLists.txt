if ( WALBERLA_BUILD_WITH_OPENMESH )
   
	waLBerla_link_geometry_to_builddir( "*.obj" )             
	waLBerla_link_files_to_builddir( "*.conf" )  

   waLBerla_add_executable( NAME ComplexGeometry FILES ComplexGeometry.cpp DEPENDS walberla::boundary walberla::core walberla::lbm walberla::mesh walberla::vtk )

	##############
	# Some tests #
	##############

   waLBerla_execute_test( NO_MODULE_LABEL NAME ComplexGeometry LABELS longrun COMMAND $<TARGET_FILE:ComplexGeometry> test.conf DEPENDS_ON_TARGETS ComplexGeometry )
endif()

