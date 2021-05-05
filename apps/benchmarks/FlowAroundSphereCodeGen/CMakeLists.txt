waLBerla_link_files_to_builddir( "*.py" )

if (WALBERLA_BUILD_WITH_CUDA)
    waLBerla_generate_target_from_python(NAME FlowAroundSphereGenerated
            FILE FlowAroundSphereCodeGen.py
            OUT_FILES FlowAroundSphereCodeGen_LbSweep.cu FlowAroundSphereCodeGen_LbSweep.h
            FlowAroundSphereCodeGen_MacroSetter.cu FlowAroundSphereCodeGen_MacroSetter.h
            FlowAroundSphereCodeGen_UBB.cu FlowAroundSphereCodeGen_UBB.h
            FlowAroundSphereCodeGen_NoSlip.cu FlowAroundSphereCodeGen_NoSlip.h
            FlowAroundSphereCodeGen_Outflow.cu FlowAroundSphereCodeGen_Outflow.h
            FlowAroundSphereCodeGen_PackInfoEven.cu FlowAroundSphereCodeGen_PackInfoEven.h
            FlowAroundSphereCodeGen_PackInfoOdd.cu FlowAroundSphereCodeGen_PackInfoOdd.h
            FlowAroundSphereCodeGen_InfoHeader.h)
    waLBerla_add_executable( NAME FlowAroundSphereCodeGen FILE FlowAroundSphereCodeGen.cpp
            DEPENDS blockforest boundary core domain_decomposition field geometry timeloop vtk FlowAroundSphereGenerated)
    set_target_properties( FlowAroundSphereCodeGen PROPERTIES CXX_VISIBILITY_PRESET hidden)
else ()
    waLBerla_generate_target_from_python(NAME FlowAroundSphereGenerated
            FILE FlowAroundSphereCodeGen.py
            OUT_FILES FlowAroundSphereCodeGen_LbSweep.cpp FlowAroundSphereCodeGen_LbSweep.h
            FlowAroundSphereCodeGen_MacroSetter.cpp FlowAroundSphereCodeGen_MacroSetter.h
            FlowAroundSphereCodeGen_UBB.cpp FlowAroundSphereCodeGen_UBB.h
            FlowAroundSphereCodeGen_NoSlip.cpp FlowAroundSphereCodeGen_NoSlip.h
            FlowAroundSphereCodeGen_Outflow.cpp FlowAroundSphereCodeGen_Outflow.h
            FlowAroundSphereCodeGen_PackInfoEven.cpp FlowAroundSphereCodeGen_PackInfoEven.h
            FlowAroundSphereCodeGen_PackInfoOdd.cpp FlowAroundSphereCodeGen_PackInfoOdd.h
            FlowAroundSphereCodeGen_InfoHeader.h)
    waLBerla_add_executable( NAME FlowAroundSphereCodeGen FILE FlowAroundSphereCodeGen.cpp
            DEPENDS blockforest boundary core domain_decomposition field geometry timeloop vtk FlowAroundSphereGenerated)
    set_target_properties( FlowAroundSphereCodeGen PROPERTIES CXX_VISIBILITY_PRESET hidden)
endif()