
waLBerla_link_files_to_builddir( "*.dat" )   
                                  
waLBerla_add_executable( NAME SchaeferTurek DEPENDS blockforest boundary core field lbm postprocessing stencil timeloop vtk sqlite )

waLBerla_execute_test( NO_MODULE_LABEL NAME SchaeferTurekTest COMMAND $<TARGET_FILE:SchaeferTurek> Test2D.dat PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )