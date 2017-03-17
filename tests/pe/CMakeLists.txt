###################################################################################################
#
# Tests for pe module
#
###################################################################################################

waLBerla_compile_test( NAME   PE_BODYFLAGS FILES BodyFlags.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_BODYFLAGS PROCESSES 8)

waLBerla_compile_test( NAME   PE_BODYITERATORS FILES BodyIterators.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_BODYITERATORS PROCESSES 2 )

waLBerla_compile_test( NAME   PE_BODYSTORAGE FILES BodyStorage.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_BODYSTORAGE )

waLBerla_compile_test( NAME   PE_CHECKVITALPARAMETERS FILES CheckVitalParameters.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_CHECKVITALPARAMETERS )

waLBerla_compile_test( NAME   PE_COLLISION FILES Collision.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_COLLISION )

waLBerla_compile_test( NAME   PE_DELETEBODY FILES DeleteBody.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_DELETEBODY_NN COMMAND $<TARGET_FILE:PE_DELETEBODY> )
waLBerla_execute_test( NAME   PE_DELETEBODY_SO COMMAND $<TARGET_FILE:PE_DELETEBODY> --syncShadowOwners )

waLBerla_compile_test( NAME   PE_DESTROYBODY FILES DestroyBody.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_DESTROYBODY )

waLBerla_compile_test( NAME   PE_FORCESYNC FILES ForceSync.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_FORCESYNC )

waLBerla_compile_test( NAME   PE_HASHGRIDS FILES HashGrids.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_HASHGRIDS_REL COMMAND $<TARGET_FILE:PE_HASHGRIDS> 10000 CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test( NAME   PE_HASHGRIDS_DBG COMMAND $<TARGET_FILE:PE_HASHGRIDS> 1000)

waLBerla_compile_test( NAME   PE_HCSITS FILES HCSITS.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_HCSITS )

waLBerla_compile_test( NAME   PE_MARSHALLING FILES Marshalling.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_MARSHALLING )

waLBerla_compile_test( NAME   PE_MATERIAL FILES Material.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_MATERIAL )

waLBerla_compile_test( NAME   PE_OVERLAP FILES Overlap.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_OVERLAP )

waLBerla_compile_test( NAME   PE_PARALLELEQUIVALENCE FILES ParallelEquivalence.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_PARALLELEQUIVALENCE PROCESSES 4 )

waLBerla_compile_test( NAME   PE_PARSEMESSAGE FILES ParseMessage.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_PARSEMESSAGE )

waLBerla_compile_test( NAME   PE_REFINEMENT FILES Refinement.cpp DEPENDS core blockforest   timeloop )
waLBerla_execute_test( NAME   PE_REFINEMENT )

waLBerla_compile_test( NAME   PE_RIGIDBODY FILES RigidBody.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_RIGIDBODY )

waLBerla_compile_test( NAME   PE_SERIALIZEDESERIALIZE FILES SerializeDeserialize.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_SERIALIZEDESERIALIZE PROCESSES 8 )

waLBerla_compile_test( NAME   PE_SHADOWCOPY FILES ShadowCopy.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   PE_SHADOWCOPY )

waLBerla_compile_test( NAME   PE_SIMPLECCD FILES SimpleCCD.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_SIMPLECCD )

waLBerla_compile_test( NAME   PE_SYNCEQUIVALENCE FILES SyncEquivalence.cpp DEPENDS core  )
#waLBerla_execute_test( NAME   PE_SYNCEQUIVALENCE COMMAND $<TARGET_FILE:PE_SYNCEQUIVALENCE> PROCESSES  8 )

waLBerla_compile_test( NAME   PE_SYNCHRONIZATION FILES Synchronization.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION01 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION03 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION09 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATION27 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATION> PROCESSES 27)

waLBerla_compile_test( NAME   PE_SYNCHRONIZATIONDELETE FILES SynchronizationDelete.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE01 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE03 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE09 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONDELETE27 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONDELETE> PROCESSES 27)

waLBerla_compile_test( NAME   PE_SYNCHRONIZATIONLARGEBODY FILES SynchronizationLargeBody.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY01 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> )
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY03 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES  3 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY09 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES  9 LABELS longrun)
waLBerla_execute_test( NAME   PE_SYNCHRONIZATIONLARGEBODY27 COMMAND $<TARGET_FILE:PE_SYNCHRONIZATIONLARGEBODY> PROCESSES 27)

waLBerla_compile_test( NAME   PE_STATICTYPEIDS FILES SetBodyTypeIDs.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_STATICTYPEIDS )

waLBerla_compile_test( NAME   PE_UNION FILES Union.cpp DEPENDS core  )
waLBerla_execute_test( NAME   PE_UNION )
