if ( WALBERLA_BUILD_WITH_OPENMESH )
        
    waLBerla_link_geometry_to_builddir( "bunny.obj" )                                
    waLBerla_add_executable( NAME MeshDistanceBenchmark DEPENDS walberla::core walberla::mesh )

	##############
	# Some tests #
	##############

	add_dependencies( waLBerlaTestsuite MeshDistanceBenchmark )
	waLBerla_execute_test( NO_MODULE_LABEL NAME MeshDistanceBenchmark COMMAND $<TARGET_FILE:MeshDistanceBenchmark> bunny.obj 10 2 )
endif()

