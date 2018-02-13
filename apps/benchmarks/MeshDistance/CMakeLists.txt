if ( WALBERLA_BUILD_WITH_OPENMESH )

	waLBerla_link_files_to_builddir( "*.obj" )                 
                                  
	waLBerla_add_executable( NAME MeshDistanceBenchmark DEPENDS core mesh )

	##############
	# Some tests #
	##############

	waLBerla_execute_test( NO_MODULE_LABEL NAME MeshDistanceBenchmark COMMAND $<TARGET_FILE:MeshDistanceBenchmark> bunny.obj 10 2 )
endif()

