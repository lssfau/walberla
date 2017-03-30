
waLBerla_link_files_to_builddir( "*.dat" )   
                                  
waLBerla_add_executable( NAME SchaeferTurek DEPENDS blockforest boundary core field lbm postprocessing stencil timeloop vtk )
