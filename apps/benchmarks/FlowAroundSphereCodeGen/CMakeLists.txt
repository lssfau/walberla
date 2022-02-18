waLBerla_link_files_to_builddir( "*.py" )

waLBerla_generate_target_from_python(NAME FlowAroundSphereGenerated
        FILE FlowAroundSphereCodeGen.py
        OUT_FILES FlowAroundSphereCodeGen_LbSweep.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_LbSweep.h
        FlowAroundSphereCodeGen_MacroSetter.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_MacroSetter.h
        FlowAroundSphereCodeGen_UBB.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_UBB.h
        FlowAroundSphereCodeGen_NoSlip.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_NoSlip.h
        FlowAroundSphereCodeGen_Outflow.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_Outflow.h
        FlowAroundSphereCodeGen_PackInfoEven.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_PackInfoEven.h
        FlowAroundSphereCodeGen_PackInfoOdd.${CODEGEN_FILE_SUFFIX} FlowAroundSphereCodeGen_PackInfoOdd.h
        FlowAroundSphereCodeGen_InfoHeader.h)

if (WALBERLA_BUILD_WITH_CUDA)
    waLBerla_add_executable( NAME FlowAroundSphereCodeGen FILE FlowAroundSphereCodeGen.cpp
            DEPENDS blockforest boundary core cuda domain_decomposition field geometry python_coupling timeloop vtk FlowAroundSphereGenerated)
else ()
    waLBerla_add_executable( NAME FlowAroundSphereCodeGen FILE FlowAroundSphereCodeGen.cpp
            DEPENDS blockforest boundary core domain_decomposition field geometry python_coupling timeloop vtk FlowAroundSphereGenerated)
endif (WALBERLA_BUILD_WITH_CUDA)