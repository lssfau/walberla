#############################################################################################################################
#
# Tests for mesh module
#
#############################################################################################################################

if ( WALBERLA_BUILD_WITH_OPENMESH )
   waLBerla_link_files_to_builddir( "*.obj" )   

   waLBerla_compile_test( FILES MatrixVectorOperationsTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MatrixVectorOperationsTest COMMAND $<TARGET_FILE:MatrixVectorOperationsTest> )
   waLBerla_compile_test( FILES MeshConversionTest.cpp DEPENDS mesh geometry )
   waLBerla_execute_test( NAME  MeshConversionTest COMMAND $<TARGET_FILE:MeshConversionTest> )

   waLBerla_compile_test( FILES MeshOperationsTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshOperationsTest COMMAND $<TARGET_FILE:MeshOperationsTest> )

   waLBerla_compile_test( FILES MeshDistanceCompareTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshDistanceCompareTest COMMAND $<TARGET_FILE:MeshDistanceCompareTest> )

   waLBerla_compile_test( FILES MeshAABBSelectionTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshAABBSelectionTestBunnyDbg COMMAND $<TARGET_FILE:MeshAABBSelectionTest> bunny.obj 2 2 )
   waLBerla_execute_test( NAME  MeshAABBSelectionTestBunny COMMAND $<TARGET_FILE:MeshAABBSelectionTest> bunny.obj 100 10 CONFIGURATIONS Release RelWithDbgInfo )
   
   waLBerla_compile_test( FILES MeshDistanceOctreeTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshDistanceOctreeTestBunny  COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> bunny.obj CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshDistanceOctreeTestW      COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> w.obj )
   waLBerla_execute_test( NAME  MeshDistanceOctreeTestCube   COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> cube.obj )
   waLBerla_execute_test( NAME  MeshDistanceOctreeTestSphere COMMAND $<TARGET_FILE:MeshDistanceOctreeTest> sphere.obj )
   
   waLBerla_compile_test( FILES MeshDistancePlausibilityTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestBunny     COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> bunny.obj  1.25   CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestBunnyLong COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> bunny.obj  0.6    CONFIGURATIONS Release RelWithDbgInfo LABELS longrun )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestW         COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> w.obj      0.005  CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestWDbg      COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> w.obj      0.1 )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestCube      COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> cube.obj   0.01 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestCubeDbg   COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> cube.obj   0.1 )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestSphere    COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> sphere.obj 0.04   CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshDistancePlausibilityTestSphereDbg COMMAND $<TARGET_FILE:MeshDistancePlausibilityTest> sphere.obj 0.1 )
   
   waLBerla_compile_test( FILES MeshContainmentOctreeTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshContainmentOctreeTestBunny  COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> bunny.obj CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshContainmentOctreeTestW      COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> w.obj )
   waLBerla_execute_test( NAME  MeshContainmentOctreeTestCube   COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> cube.obj )
   waLBerla_execute_test( NAME  MeshContainmentOctreeTestSphere COMMAND $<TARGET_FILE:MeshContainmentOctreeTest> sphere.obj )
   
   waLBerla_compile_test( FILES MeshVTKTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshVTKTest  COMMAND $<TARGET_FILE:MeshVTKTest> sphere.obj )
   
   waLBerla_compile_test( FILES DistributedMeshVTKTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  DistributedMeshVTKTest  COMMAND $<TARGET_FILE:DistributedMeshVTKTest> sphere.obj )
   
   waLBerla_compile_test( FILES PeVTKMeshWriterTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  PeVTKMeshWriterTest  COMMAND $<TARGET_FILE:PeVTKMeshWriterTest> )
   
   waLBerla_compile_test( FILES MeshPeRaytracing.cpp DEPENDS mesh pe )
   waLBerla_execute_test( NAME  MeshPeRaytracing )
   
   waLBerla_compile_test( FILES MeshAABBIntersectionTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshAABBIntersectionTestDbg COMMAND $<TARGET_FILE:MeshAABBIntersectionTest> 10                                          )
   waLBerla_execute_test( NAME  MeshAABBIntersectionTest    COMMAND $<TARGET_FILE:MeshAABBIntersectionTest> 10000 CONFIGURATIONS Release RelWithDbgInfo )
   
   waLBerla_compile_test( FILES MeshInitilizationTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshInitializationTestDbg1P1 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 1 1  LABELS longrun )
   waLBerla_execute_test( NAME  MeshInitializationTestDbg2P1 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 2 12 LABELS longrun )
   waLBerla_execute_test( NAME  MeshInitializationTestDbg2P2 COMMAND $<TARGET_FILE:MeshInitilizationTest> cube.obj 2 12 PROCESSES 2    )
   
   waLBerla_execute_test( NAME  MeshInitializationTest1P    COMMAND $<TARGET_FILE:MeshInitilizationTest> bunny.obj 6 40 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshInitializationTest4P    COMMAND $<TARGET_FILE:MeshInitilizationTest> bunny.obj 6 40 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )

   waLBerla_compile_test( FILES MeshInitilizationHelperTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg1P1 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 1 1              )
   waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg2P1 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 2 12             )
   waLBerla_execute_test( NAME  MeshInitializationHelperTestDbg2P2 COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> cube.obj 2 12 PROCESSES 2 )
   
   waLBerla_execute_test( NAME  MeshInitializationHelperTest1P    COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> bunny.obj 6 40 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshInitializationHelperTest4P    COMMAND $<TARGET_FILE:MeshInitilizationHelperTest> bunny.obj 6 40 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
   
   waLBerla_compile_test( FILES MeshBlockExclusionTest.cpp DEPENDS mesh )
   waLBerla_execute_test( NAME  MeshBlockExclusionTestBunny  COMMAND $<TARGET_FILE:MeshBlockExclusionTest> bunny.obj  10000 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshBlockExclusionTestW      COMMAND $<TARGET_FILE:MeshBlockExclusionTest> w.obj      1000 )
   waLBerla_execute_test( NAME  MeshBlockExclusionTestCube   COMMAND $<TARGET_FILE:MeshBlockExclusionTest> cube.obj   1000 )
   waLBerla_execute_test( NAME  MeshBlockExclusionTestSphere COMMAND $<TARGET_FILE:MeshBlockExclusionTest> sphere.obj 1000 )
   
   waLBerla_execute_test( NAME  MeshBlockExclusionTest2PBunny  COMMAND $<TARGET_FILE:MeshBlockExclusionTest> bunny.obj  10000 PROCESSES 2 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  MeshBlockExclusionTest2PW      COMMAND $<TARGET_FILE:MeshBlockExclusionTest> w.obj      1000 PROCESSES 2 )
   waLBerla_execute_test( NAME  MeshBlockExclusionTest2PCube   COMMAND $<TARGET_FILE:MeshBlockExclusionTest> cube.obj   1000 PROCESSES 2 )
   waLBerla_execute_test( NAME  MeshBlockExclusionTest2PSphere COMMAND $<TARGET_FILE:MeshBlockExclusionTest> sphere.obj 1000 PROCESSES 2 )

   waLBerla_compile_test( FILES OpenMeshBufferTypeExtensionsTest.cpp DEPENDS core mesh )
   waLBerla_execute_test( NAME  OpenMeshBufferTypeExtensionsTest COMMAND $<TARGET_FILE:OpenMeshBufferTypeExtensionsTest> )

   waLBerla_compile_test( FILES QHullTest.cpp DEPENDS core mesh )
   waLBerla_execute_test( NAME  QHullTest             COMMAND $<TARGET_FILE:QHullTest>             10000 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  QHullTestFloatMesh    COMMAND $<TARGET_FILE:QHullTest> --floatMesh 10000 CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME  QHullTestDbg          COMMAND $<TARGET_FILE:QHullTest>             10 )
   waLBerla_execute_test( NAME  QHullTestFloatMeshDbg COMMAND $<TARGET_FILE:QHullTest> --floatMesh 10 )
   
   waLBerla_compile_test( FILES NumericIntegrationTest.cpp DEPENDS core mesh )
   waLBerla_execute_test( NAME  NumericIntegrationTestCube   COMMAND $<TARGET_FILE:NumericIntegrationTest> cube.obj )
   waLBerla_execute_test( NAME  NumericIntegrationTestSphere COMMAND $<TARGET_FILE:NumericIntegrationTest> sphere.obj )
   waLBerla_execute_test( NAME  NumericIntegrationTestBunny  COMMAND $<TARGET_FILE:NumericIntegrationTest> bunny.obj CONFIGURATIONS Release RelWithDbgInfo )

   waLBerla_compile_test( FILES MeshMarshalling.cpp DEPENDS core mesh )
   waLBerla_execute_test( NAME  MeshMarshalling )
endif()
