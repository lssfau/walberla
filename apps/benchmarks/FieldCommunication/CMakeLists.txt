
waLBerla_link_files_to_builddir( "*.dat" )
waLBerla_link_files_to_builddir( "*.py" )


waLBerla_add_executable ( NAME FieldCommunication
      DEPENDS walberla::blockforest walberla::core walberla::domain_decomposition walberla::field walberla::postprocessing walberla::sqlite walberla::python_coupling )
