###################################################################################################
#
# Tests for walberla::pe module
#
###################################################################################################

waLBerla_link_files_to_builddir( *.cfg )
waLBerla_link_files_to_builddir( *.sbf )

waLBerla_compile_test( NAME PE_BODYFLAGS FILES BodyFlags.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_BODYFLAGS PROCESSES 8)

waLBerla_compile_test( NAME PE_BODYITERATORS FILES BodyIterators.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_BODYITERATORS PROCESSES 2 )

waLBerla_compile_test( NAME PE_BODYSTORAGE FILES BodyStorage.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_BODYSTORAGE )

waLBerla_compile_test( NAME PE_CALLBACK FILES Callback.cpp DEPENDS walberla::blockforest walberla::core walberla::domain_decomposition )
waLBerla_execute_test( NAME   PE_CALLBACK PROCESSES 2 )

waLBerla_compile_test( NAME PE_CHECKVITALPARAMETERS FILES CheckVitalParameters.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_CHECKVITALPARAMETERS )

waLBerla_compile_test( NAME PE_COLLISION FILES Collision.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_COLLISION )

waLBerla_compile_test( NAME PE_COLLISIONTOBIASGJK FILES CollisionTobiasGJK.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_COLLISIONTOBIASGJK )

waLBerla_compile_test( NAME PE_CREATEWORLD FILES CreateWorld.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_CREATEWORLD )

waLBerla_compile_test( NAME PE_DELETEBODY FILES DeleteBody.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_DELETEBODY_NN COMMAND $<TARGET_FILE:PE_DELETEBODY> )
waLBerla_execute_test( NAME   PE_DELETEBODY_SO COMMAND $<TARGET_FILE:PE_DELETEBODY> --syncShadowOwners )

waLBerla_compile_test( NAME PE_DESTROYBODY FILES DestroyBody.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_DESTROYBODY )

waLBerla_compile_test( NAME PE_DOCUMENTATIONSNIPPETS FILES PeDocumentationSnippets.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_DOCUMENTATIONSNIPPETS )

waLBerla_compile_test( NAME PE_DYNAMICREFINEMENT FILES DynamicRefinement.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_DYNAMICREFINEMENT )

waLBerla_compile_test( NAME PE_FORCESYNC FILES ForceSync.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_FORCESYNC )

waLBerla_compile_test( NAME PE_HASHGRIDS FILES HashGrids.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_HASHGRIDS_REL COMMAND $<TARGET_FILE:PE_HASHGRIDS> 10000 CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test( NAME   PE_HASHGRIDS_DBG COMMAND $<TARGET_FILE:PE_HASHGRIDS> 1000)

waLBerla_compile_test( NAME PE_HCSITS FILES HCSITS.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_HCSITS )

waLBerla_compile_test( NAME PE_LOADFROMCONFIG FILES LoadFromConfig.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_LOADFROMCONFIG COMMAND $<TARGET_FILE:PE_LOADFROMCONFIG> ${CMAKE_CURRENT_SOURCE_DIR}/LoadFromConfig.cfg )

waLBerla_compile_test( NAME PE_MARSHALLING FILES Marshalling.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_MARSHALLING )

waLBerla_compile_test( NAME PE_MATERIAL FILES Material.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_MATERIAL )

waLBerla_compile_test( NAME PE_MINMAXREFINEMENT FILES MinMaxRefinement.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_MINMAXREFINEMENT PROCESSES 8 )

waLBerla_compile_test( NAME PE_OVERLAP FILES Overlap.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_OVERLAP )

waLBerla_compile_test( NAME PE_PARALLELEQUIVALENCE FILES ParallelEquivalence.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_PARALLELEQUIVALENCE PROCESSES 4 )

if( WALBERLA_BUILD_WITH_PARMETIS )
   waLBerla_compile_test( NAME PE_PARMETIS FILES ParMetis.cpp DEPENDS walberla::core walberla::blockforest )
   waLBerla_execute_test( NAME   PE_PARMETIS PROCESSES 64 )
endif()

waLBerla_compile_test( NAME PE_PARSEMESSAGE FILES ParseMessage.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_PARSEMESSAGE )

waLBerla_compile_test( NAME PE_REFINEMENT FILES Refinement.cpp DEPENDS walberla::core walberla::blockforest walberla::timeloop )
waLBerla_execute_test( NAME   PE_REFINEMENT )

waLBerla_compile_test( NAME PE_REFINEMENTWITHSPAREPROCESSES FILES RefinementWithSpareProcesses.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_REFINEMENTWITHSPAREPROCESSES PROCESSES 27 )

waLBerla_compile_test( NAME PE_RIGIDBODY FILES RigidBody.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_RIGIDBODY )

waLBerla_compile_test( NAME PE_SERIALIZEDESERIALIZE FILES SerializeDeserialize.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   PE_SERIALIZEDESERIALIZE01 COMMAND $<TARGET_FILE:PE_SERIALIZEDESERIALIZE> PROCESSES 1 )
waLBerla_execute_test( NAME   PE_SERIALIZEDESERIALIZE02 COMMAND $<TARGET_FILE:PE_SERIALIZEDESERIALIZE> PROCESSES 2 )
waLBerla_execute_test( NAME   PE_SERIALIZEDESERIALIZE04 COMMAND $<TARGET_FILE:PE_SERIALIZEDESERIALIZE> PROCESSES 4 )
waLBerla_execute_test( NAME   PE_SERIALIZEDESERIALIZE08 COMMAND $<TARGET_FILE:PE_SERIALIZEDESERIALIZE> PROCESSES 8 )
if( WALBERLA_BUILD_WITH_MPI )
set_property( TEST PE_SERIALIZEDESERIALIZE02 PROPERTY DEPENDS PE_SERIALIZEDESERIALIZE01 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST PE_SERIALIZEDESERIALIZE04 PROPERTY DEPENDS PE_SERIALIZEDESERIALIZE02 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST PE_SERIALIZEDESERIALIZE08 PROPERTY DEPENDS PE_SERIALIZEDESERIALIZE04 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
endif()

waLBerla_compile_test( NAME PE_SHADOWCOPY FILES ShadowCopy.cpp DEPENDS walberla::core walberla::blockforest walberla::domain_decomposition )
waLBerla_execute_test( NAME   PE_SHADOWCOPY_NN COMMAND $<TARGET_FILE:PE_SHADOWCOPY> )
waLBerla_execute_test( NAME   PE_SHADOWCOPY_SO COMMAND $<TARGET_FILE:PE_SHADOWCOPY> --syncShadowOwners )

waLBerla_compile_test( NAME PE_SIMPLECCD FILES SimpleCCD.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_SIMPLECCD )

waLBerla_compile_test( NAME PE_SYNCEQUIVALENCE FILES SyncEquivalence.cpp DEPENDS walberla::core )
#waLBerla_execute_test( NAME   PE_SYNCEQUIVALENCE COMMAND $<TARGET_FILE:PE_SYNCEQUIVALENCE> PROCESSES  8 )

waLBerla_compile_test( NAME PE_SYNCHRONIZATION FILES Synchronization.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION01 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION03 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION09 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION27 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES 27)

waLBerla_compile_test( NAME PE_SYNCHRONIZATIONDELETE FILES SynchronizationDelete.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE01_NN COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE03_NN COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE09_NN COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE27_NN COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES 27)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE01_SO COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> --syncShadowOwners )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE03_SO COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> --syncShadowOwners PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE09_SO COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> --syncShadowOwners PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE27_SO COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> --syncShadowOwners PROCESSES 27)

waLBerla_compile_test( NAME PE_SYNCHRONIZATIONLARGEBODY FILES SynchronizationLargeBody.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY01 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY03 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY09 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY27 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES 27)

waLBerla_compile_test( NAME PE_STATICTYPEIDS FILES SetBodyTypeIDs.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_STATICTYPEIDS )

waLBerla_compile_test( NAME PE_UNION FILES Union.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_UNION )

waLBerla_compile_test( NAME PE_UNIONBEHAVIOR FILES UnionBehavior.cpp DEPENDS walberla::core )
if( WALBERLA_DOUBLE_ACCURACY )
waLBerla_execute_test( NAME   PE_UNIONBEHAVIOR )
endif()

waLBerla_compile_test( NAME PE_RAYTRACING FILES Raytracing.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_RAYTRACING )

waLBerla_compile_test( NAME PE_RAYTRACING_INTERSECTION FILES RaytracingIntersection.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_RAYTRACING_INTERSECTION )

waLBerla_compile_test( NAME PE_RAYTRACING_SPHERE FILES RaytracingSphere.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_RAYTRACING_SPHERE )

waLBerla_compile_test( NAME PE_RAYTRACING_HASHGRIDS FILES RaytracingHashGrids.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_RAYTRACING_HASHGRIDS )

waLBerla_compile_test( NAME PE_VOLUMEINERTIA FILES VolumeInertia.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME   PE_VOLUMEINERTIA CONFIGURATIONS Release RelWithDbgInfo)
