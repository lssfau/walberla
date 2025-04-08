waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )

waLBerla_generate_target_from_python(NAME NonUniformGridCPUGenerated
        FILE NonUniformGridCPU.py
        OUT_FILES NonUniformGridCPUStorageSpecification.h NonUniformGridCPUStorageSpecification.cpp
        NonUniformGridCPUSweepCollection.h NonUniformGridCPUSweepCollection.cpp
        NoSlip.h NoSlip.cpp
        UBB.h UBB.cpp
        NonUniformGridCPUBoundaryCollection.h
        NonUniformGridCPUInfoHeader.h)

waLBerla_add_executable( NAME NonUniformGridCPU
                         FILES NonUniformGridCPU.cpp LdcSetup.h GridGeneration.h
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::geometry walberla::lbm_generated walberla::python_coupling walberla::timeloop walberla::vtk NonUniformGridCPUGenerated )
