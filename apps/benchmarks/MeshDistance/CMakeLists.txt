if ( WALBERLA_BUILD_WITH_OPENMESH )
        
    waLBerla_link_geometry_to_builddir( "bunny.obj" )                                
    add_executable ( MeshDistanceBenchmark MeshDistanceBenchmark.cpp )
    target_link_libraries ( MeshDistanceBenchmark walberla::core walberla::mesh )

	##############
	# Some tests #
	##############

	add_dependencies( waLBerlaTestsuite MeshDistanceBenchmark )
	waLBerla_execute_test( NO_MODULE_LABEL NAME MeshDistanceBenchmark COMMAND $<TARGET_FILE:MeshDistanceBenchmark> bunny.obj 10 2 )
endif()

