waLBerla_link_files_to_builddir( *.prm )

waLBerla_generate_target_from_python(NAME MeshRefinementExampleCodeGen
    FILE MeshRefinementExample.py
    OUT_FILES NoSlipBC.h NoSlipBC.cpp
    UBBBC.h UBBBC.cpp
    OutflowBC.h OutflowBC.cpp
    FreeSlipBC.h FreeSlipBC.cpp
    MeshRefinementExampleStorageSpecification.h MeshRefinementExampleStorageSpecification.cpp
    MeshRefinementExampleSweepCollection.h MeshRefinementExampleSweepCollection.cpp
    MeshRefinementExampleBoundaryCollection.h
    InfoHeader.h
)

add_executable ( MeshRefinementExample MeshRefinementExample.cpp)
target_link_libraries ( MeshRefinementExample  walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::geometry walberla::lbm_generated
        walberla::timeloop walberla::vtk MeshRefinementExampleCodeGen )