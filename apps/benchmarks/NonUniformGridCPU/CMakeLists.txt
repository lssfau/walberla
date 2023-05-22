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
                         FILES NonUniformGridCPU.cpp
                         DEPENDS blockforest boundary core domain_decomposition field geometry python_coupling timeloop vtk NonUniformGridCPUGenerated )