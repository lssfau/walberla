
waLBerla_link_files_to_builddir( "*.dat" )
waLBerla_link_files_to_builddir( "*.py" )


waLBerla_add_executable ( NAME FieldCommunication
                          DEPENDS blockforest core domain_decomposition field postprocessing sqlite python_coupling )
