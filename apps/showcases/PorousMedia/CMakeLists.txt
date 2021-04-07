waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.txt )

waLBerla_add_executable (
        NAME PackedBedCreation
        FILES PackedBedCreation.cpp
        DEPENDS blockforest core field geometry timeloop pe pe_coupling vtk)

waLBerla_generate_target_from_python(
        NAME GeneratedLatticeModel
        FILE PorousMediaCumulantLBMKernel.py
        OUT_FILES LbCodeGen_LatticeModel.cpp LbCodeGen_LatticeModel.h)

waLBerla_add_executable (
        NAME PorousMedia
        FILES PorousMediaCumulantLBMKernel.py PorousMedia.cpp
        DEPENDS blockforest core field lbm geometry timeloop pe pe_coupling GeneratedLatticeModel )

