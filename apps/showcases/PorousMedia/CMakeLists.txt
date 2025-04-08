waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.txt )

waLBerla_add_executable (
        NAME PackedBedCreation
        FILES PackedBedCreation.cpp
      DEPENDS walberla::blockforest walberla::core walberla::field walberla::geometry walberla::timeloop walberla::pe walberla::pe_coupling walberla::vtk )

waLBerla_generate_target_from_python(
        NAME GeneratedLatticeModel
        FILE PorousMediaCumulantLBMKernel.py
        OUT_FILES LbCodeGen_LatticeModel.cpp LbCodeGen_LatticeModel.h)

waLBerla_add_executable (
        NAME PorousMedia
        FILES PorousMediaCumulantLBMKernel.py PorousMedia.cpp
      DEPENDS walberla::blockforest walberla::core walberla::field walberla::lbm walberla::geometry walberla::timeloop walberla::pe walberla::pe_coupling GeneratedLatticeModel )

