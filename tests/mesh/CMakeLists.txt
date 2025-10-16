#############################################################################################################################
#
# Tests for walberla::mesh module
#
#############################################################################################################################

if ( WALBERLA_BUILD_WITH_OPENMESH )
    waLBerla_link_files_to_builddir( "w.obj" )
    waLBerla_link_geometry_to_builddir( "*.obj" )

    waLBerla_add_test_executable ( MatrixVectorOperationsTest MatrixVectorOperationsTest.cpp )
    target_link_libraries ( MatrixVectorOperationsTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MatrixVectorOperationsTest COMMAND $<TARGET_FILE:MatrixVectorOperationsTest> )

    waLBerla_add_test_executable ( MeshConversionTest MeshConversionTest.cpp )
    target_link_libraries ( MeshConversionTest walberla::mesh walberla::geometry walberla::mesh )
    waLBerla_execute_test( NAME  MeshConversionTest COMMAND $<TARGET_FILE:MeshConversionTest> )

    waLBerla_add_test_executable ( MeshOperationsTest MeshOperationsTest.cpp )
    target_link_libraries ( MeshOperationsTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshOperationsTest COMMAND $<TARGET_FILE:MeshOperationsTest> )

    waLBerla_add_test_executable ( MeshDistanceCompareTest MeshDistanceCompareTest.cpp )
    target_link_libraries ( MeshDistanceCompareTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshDistanceCompareTest COMMAND $<TARGET_FILE:MeshDistanceCompareTest> )

    waLBerla_add_test_executable ( MeshAABBSelectionTest MeshAABBSelectionTest.cpp )
    target_link_libraries ( MeshAABBSelectionTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshAABBSelectionTestBunnyDbg COMMAND $<TARGET_FILE:MeshAABBSelectionTest> bunny.obj 2 2 )
    waLBerla_execute_test( NAME  MeshAABBSelectionTestBunny COMMAND $<TARGET_FILE:MeshAABBSelectionTest> bunny.obj 100 10 )

    waLBerla_add_test_executable ( MeshDistanceOctreeTest MeshDistanceOctreeTest.cpp )
    target_link_libraries ( MeshDistanceOctreeTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshDistanceOctreeTestBunny  COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> bunny.obj )
    waLBerla_execute_test( NAME  MeshDistanceOctreeTestW      COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> w.obj )
    waLBerla_execute_test( NAME  MeshDistanceOctreeTestCube   COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> cube.obj )
    waLBerla_execute_test( NAME  MeshDistanceOctreeTestSphere COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> sphere.obj )

    waLBerla_add_test_executable ( MeshDistanceMollerTrumboreTest MeshDistanceMollerTrumboreTest.cpp )
    target_link_libraries ( MeshDistanceMollerTrumboreTest walberla::mesh walberla::mesh )

    waLBerla_add_test_executable ( MeshDistanceCallbackMollerTrumboreTest MeshDistanceCallbackMollerTrumboreTest.cpp )
    target_link_libraries ( MeshDistanceCallbackMollerTrumboreTest walberla::mesh walberla::mesh )

    if( NOT WALBERLA_OPTIMIZE_FOR_LOCALHOST )
        waLBerla_execute_test( NAME  MeshDistanceMollerTrumboreTestCube COMMAND $<TARGET_FILE:MeshDistanceMollerTrumboreTest>  cube.obj )

        waLBerla_execute_test( NAME  MeshDistanceCallbackMollerTrumboreTestCube COMMAND $<TARGET_FILE:MeshDistanceCallbackMollerTrumboreTest>  cube.obj )
        if ( WALBERLA_DOUBLE_ACCURACY )
            waLBerla_execute_test(
            NAME  MeshDistanceCallbackMollerTrumboreTestSphere
            COMMAND $<TARGET_FILE:MeshDistanceCallbackMollerTrumboreTest>  sphere.obj
            LABELS longrun
            )
        endif ()
    endif ()

    waLBerla_add_test_executable ( MeshDistancePlausibilityTest MeshDistancePlausibilityTest.cpp )
    target_link_libraries ( MeshDistancePlausibilityTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestBunny     COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> bunny.obj  1.25   )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestBunnyLong COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> bunny.obj  0.6    LABELS longrun )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestW         COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> w.obj      0.005  )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestWDbg      COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> w.obj      0.1 )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestCube      COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> cube.obj   0.01 )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestCubeDbg   COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> cube.obj   0.1 )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestSphere    COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> sphere.obj 0.04   )
    waLBerla_execute_test( NAME  MeshDistancePlausibilityTestSphereDbg COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> sphere.obj 0.1 )

    waLBerla_add_test_executable ( MeshContainmentOctreeTest MeshContainmentOctreeTest.cpp )
    target_link_libraries ( MeshContainmentOctreeTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshContainmentOctreeTestBunny  COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> bunny.obj )
    waLBerla_execute_test( NAME  MeshContainmentOctreeTestW      COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> w.obj )
    waLBerla_execute_test( NAME  MeshContainmentOctreeTestCube   COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> cube.obj )
    waLBerla_execute_test( NAME  MeshContainmentOctreeTestSphere COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> sphere.obj )

    waLBerla_add_test_executable ( MeshVTKTest MeshVTKTest.cpp )
    target_link_libraries ( MeshVTKTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshVTKTest  COMMAND $<TARGET_FILE:MeshVTKTest> sphere.obj )

    waLBerla_add_test_executable ( DistributedMeshVTKTest DistributedMeshVTKTest.cpp )
    target_link_libraries ( DistributedMeshVTKTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  DistributedMeshVTKTest  COMMAND $<TARGET_FILE:DistributedMeshVTKTest> sphere.obj )

    waLBerla_add_test_executable ( MeshAABBIntersectionTest MeshAABBIntersectionTest.cpp )
    target_link_libraries ( MeshAABBIntersectionTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshAABBIntersectionTestDbg COMMAND $<TARGET_FILE:MeshAABBIntersectionTest> 10                                          )
    waLBerla_execute_test( NAME  MeshAABBIntersectionTest    COMMAND $<TARGET_FILE:MeshAABBIntersectionTest> 10000 )

    waLBerla_add_test_executable ( MeshInitilizationTest MeshInitilizationTest.cpp )
    target_link_libraries ( MeshInitilizationTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshInitializationTestDbg1P1 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 1 1  LABELS longrun )
    waLBerla_execute_test( NAME  MeshInitializationTestDbg2P1 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 2 12 LABELS longrun )
    waLBerla_execute_test( NAME  MeshInitializationTestDbg2P2 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 2 12 PROCESSES 2 )

    waLBerla_execute_test( NAME  MeshInitializationTest1P    COMMAND $<TARGET_FILE:MeshInitilizationTest> bunny.obj 6 40 LABELS longrun )
    waLBerla_execute_test( NAME  MeshInitializationTest4P    COMMAND $<TARGET_FILE:MeshInitilizationTest> bunny.obj 6 40 PROCESSES 4 LABELS longrun )

    waLBerla_add_test_executable ( MeshInitilizationHelperTest MeshInitilizationHelperTest.cpp )
    target_link_libraries ( MeshInitilizationHelperTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg1P1 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 1 1              )
    waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg2P1 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 2 12             )
    waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg2P2 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 2 12 PROCESSES 2 )

    waLBerla_execute_test( NAME  MeshInitializationHelperTest1P    COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> bunny.obj 6 40 )
    waLBerla_execute_test( NAME  MeshInitializationHelperTest4P    COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> bunny.obj 6 40 PROCESSES 4 )

    waLBerla_add_test_executable ( MeshBlockExclusionTest MeshBlockExclusionTest.cpp )
    target_link_libraries ( MeshBlockExclusionTest walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  MeshBlockExclusionTestBunny  COMMAND $<TARGET_FILE:MeshBlockExclusionTest> bunny.obj  10000 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTestW      COMMAND $<TARGET_FILE:MeshBlockExclusionTest> w.obj      1000 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTestCube   COMMAND $<TARGET_FILE:MeshBlockExclusionTest> cube.obj   1000 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTestSphere COMMAND $<TARGET_FILE:MeshBlockExclusionTest> sphere.obj 1000 )

    waLBerla_execute_test( NAME  MeshBlockExclusionTest2PBunny  COMMAND $<TARGET_FILE:MeshBlockExclusionTest> bunny.obj  10000 PROCESSES 2 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTest2PW      COMMAND $<TARGET_FILE:MeshBlockExclusionTest> w.obj      1000 PROCESSES 2 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTest2PCube   COMMAND $<TARGET_FILE:MeshBlockExclusionTest> cube.obj   1000 PROCESSES 2 )
    waLBerla_execute_test( NAME  MeshBlockExclusionTest2PSphere COMMAND $<TARGET_FILE:MeshBlockExclusionTest> sphere.obj 1000 PROCESSES 2 )

    waLBerla_add_test_executable ( OpenMeshBufferTypeExtensionsTest OpenMeshBufferTypeExtensionsTest.cpp )
    target_link_libraries ( OpenMeshBufferTypeExtensionsTest walberla::core walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  OpenMeshBufferTypeExtensionsTest COMMAND $<TARGET_FILE:OpenMeshBufferTypeExtensionsTest> )

    waLBerla_add_test_executable ( QHullTest QHullTest.cpp )
    target_link_libraries ( QHullTest walberla::core walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  QHullTest             COMMAND $<TARGET_FILE:QHullTest>             10000 )
    waLBerla_execute_test( NAME  QHullTestFloatMesh    COMMAND $<TARGET_FILE:QHullTest> --floatMesh 10000 )
    waLBerla_execute_test( NAME  QHullTestDbg          COMMAND $<TARGET_FILE:QHullTest>             10 )
    waLBerla_execute_test( NAME  QHullTestFloatMeshDbg COMMAND $<TARGET_FILE:QHullTest> --floatMesh 10 )

    waLBerla_add_test_executable ( NumericIntegrationTest NumericIntegrationTest.cpp )
    target_link_libraries ( NumericIntegrationTest walberla::core walberla::mesh walberla::mesh )
    waLBerla_execute_test( NAME  NumericIntegrationTestCube   COMMAND $<TARGET_FILE:NumericIntegrationTest> cube.obj )
    waLBerla_execute_test( NAME  NumericIntegrationTestSphere COMMAND $<TARGET_FILE:NumericIntegrationTest> sphere.obj )
    waLBerla_execute_test( NAME  NumericIntegrationTestBunny  COMMAND $<TARGET_FILE:NumericIntegrationTest> bunny.obj )
endif()
