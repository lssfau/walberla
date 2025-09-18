
waLBerla_link_files_to_builddir( "*.dat" )
waLBerla_link_files_to_builddir( "*.py" )


add_executable ( FieldCommunication FieldCommunication.cpp )
target_link_libraries ( FieldCommunication walberla::blockforest walberla::core walberla::domain_decomposition walberla::field walberla::postprocessing walberla::sqlite walberla::python_coupling )
