waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.emmy )

waLBerla_add_executable( NAME HourGlass 
                         FILES HourGlass.cpp  
                         DEPENDS blockforest domain_decomposition pe postprocessing vtk  )
