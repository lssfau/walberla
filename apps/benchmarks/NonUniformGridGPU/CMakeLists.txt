waLBerla_link_files_to_builddir( "*.prm" )
waLBerla_link_files_to_builddir( "*.py" )
waLBerla_link_files_to_builddir( "simulation_setup" )

waLBerla_generate_target_from_python(NAME NonUniformGridGPUGenerated
        FILE NonUniformGridGPU.py
        OUT_FILES NonUniformGridGPUStorageSpecification.h NonUniformGridGPUStorageSpecification.${CODEGEN_FILE_SUFFIX}
        NonUniformGridGPUSweepCollection.h NonUniformGridGPUSweepCollection.${CODEGEN_FILE_SUFFIX}
        NoSlip.h NoSlip.${CODEGEN_FILE_SUFFIX}
        UBB.h UBB.${CODEGEN_FILE_SUFFIX}
        NonUniformGridGPUBoundaryCollection.h
        NonUniformGridGPUInfoHeader.h)
waLBerla_add_executable( NAME NonUniformGridGPU
                         FILES NonUniformGridGPU.cpp LdcSetup.h
                         DEPENDS blockforest boundary core gpu domain_decomposition field geometry lbm_generated python_coupling timeloop vtk NonUniformGridGPUGenerated )