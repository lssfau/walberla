###################################################################################################
#
# Tests for MESA_PD
#
###################################################################################################

waLBerla_link_files_to_builddir(*.prm)

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_AnalyticCollisionFunctions collision_detection/AnalyticCollisionFunctions.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_AnalyticCollisionFunctions walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticCollisionFunctions )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_AnalyticContactDetection collision_detection/AnalyticContactDetection.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_AnalyticContactDetection walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticContactDetection )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_BoxSupport collision_detection/BoxSupport.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_BoxSupport walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_BoxSupport )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_EllipsoidSupport collision_detection/EllipsoidSupport.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_EllipsoidSupport walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_EllipsoidSupport )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_EPA collision_detection/EPA.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_EPA walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_EPA )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_GeneralContactDetection collision_detection/GeneralContactDetection.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_GeneralContactDetection walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GeneralContactDetection )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_GJK collision_detection/GJK.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_GJK walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GJK )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_GJKEPA collision_detection/GJK_EPA.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_GJKEPA walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GJKEPA )

waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_SphereSupport collision_detection/SphereSupport.cpp )
target_link_libraries ( MESA_PD_COLLISIONDETECTION_SphereSupport walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_SphereSupport )

waLBerla_add_test_executable ( MESA_PD_COMMON_IntersectionRatio common/IntersectionRatio.cpp )
target_link_libraries ( MESA_PD_COMMON_IntersectionRatio walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COMMON_IntersectionRatio )

waLBerla_add_test_executable ( MESA_PD_COMMON_ContainsPoint common/ContainsPoint.cpp )
target_link_libraries ( MESA_PD_COMMON_ContainsPoint walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_COMMON_ContainsPoint )

waLBerla_add_test_executable ( MESA_PD_ContactDetection ContactDetection.cpp )
target_link_libraries ( MESA_PD_ContactDetection walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_ContactDetection PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_Data_ContactHistory data/ContactHistory.cpp )
target_link_libraries ( MESA_PD_Data_ContactHistory walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_ContactHistory )

waLBerla_add_test_executable ( MESA_PD_Data_Flags data/Flags.cpp )
target_link_libraries ( MESA_PD_Data_Flags walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_Flags )

waLBerla_add_test_executable ( MESA_PD_Data_HashGridsVsBruteForce data/HashGridsVsBruteForce.cpp )
target_link_libraries ( MESA_PD_Data_HashGridsVsBruteForce walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_HashGridsVsBruteForce PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Data_LinkedCells data/LinkedCells.cpp )
target_link_libraries ( MESA_PD_Data_LinkedCells walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_LinkedCells )

waLBerla_add_test_executable ( MESA_PD_Data_ParticleStorage data/ParticleStorage.cpp )
target_link_libraries ( MESA_PD_Data_ParticleStorage walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_ParticleStorage )

waLBerla_add_test_executable ( MESA_PD_Data_SparseLinkedCells data/SparseLinkedCells.cpp )
target_link_libraries ( MESA_PD_Data_SparseLinkedCells walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Data_SparseLinkedCells )

waLBerla_add_test_executable ( MESA_PD_Domain_BlockForestDomain domain/BlockForestDomain.cpp )
target_link_libraries ( MESA_PD_Domain_BlockForestDomain walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestDomain )

waLBerla_add_test_executable ( MESA_PD_Domain_BlockForestSync domain/BlockForestSync.cpp )
target_link_libraries ( MESA_PD_Domain_BlockForestSync walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync2 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 2 )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync4 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 4 )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync8 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_Domain_BlockForestSyncPeriodic domain/BlockForestSyncPeriodic.cpp )
target_link_libraries ( MESA_PD_Domain_BlockForestSyncPeriodic walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSyncPeriodic PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_Domain_DistanceCalculation domain/DistanceCalculation.cpp )
target_link_libraries ( MESA_PD_Domain_DistanceCalculation walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_DistanceCalculation )

waLBerla_add_test_executable ( MESA_PD_Domain_DynamicRefinement domain/DynamicRefinement.cpp )
target_link_libraries ( MESA_PD_Domain_DynamicRefinement walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_DynamicRefinement PROCESSES 8)

waLBerla_add_test_executable ( MESA_PD_Domain_InfiniteDomain domain/InfiniteDomain.cpp )
target_link_libraries ( MESA_PD_Domain_InfiniteDomain walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_InfiniteDomain )

waLBerla_add_test_executable ( MESA_PD_Domain_SerializeDeserialize domain/SerializeDeserialize.cpp )
target_link_libraries ( MESA_PD_Domain_SerializeDeserialize walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Domain_SerializeDeserialize PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_DropTestAnalytic DropTestAnalytic.cpp )
target_link_libraries ( MESA_PD_DropTestAnalytic walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_DropTestAnalytic )

waLBerla_add_test_executable ( MESA_PD_DropTestGeneral DropTestGeneral.cpp )
target_link_libraries ( MESA_PD_DropTestGeneral walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_DropTestGeneral )

waLBerla_add_test_executable ( MESA_PD_Kernel_ClearNextNeighborSync kernel/ClearNextNeighborSync.cpp )
target_link_libraries ( MESA_PD_Kernel_ClearNextNeighborSync walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ClearNextNeighborSync PROCESSES 2 )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_AnisotropicVDWContact kernel/cnt/AnisotropicVDWContact.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_AnisotropicVDWContact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_AnisotropicVDWContact )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_IntegratedVDWContact kernel/cnt/IntegratedVDWContact.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_IntegratedVDWContact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_IntegratedVDWContact )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_IsotropicVDWContact kernel/cnt/IsotropicVDWContact.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_IsotropicVDWContact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_IsotropicVDWContact )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_VBondContact kernel/cnt/VBondContact.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_VBondContact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_VBondContact )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_VBondContactIntegration kernel/cnt/VBondContactIntegration.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_VBondContactIntegration walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME MESA_PD_Kernel_CNT_VBondContactIntegration )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_ViscousDamping kernel/cnt/ViscousDamping.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_ViscousDamping walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_ViscousDamping )

waLBerla_add_test_executable ( MESA_PD_Kernel_CNT_WallContact kernel/cnt/WallContact.test.cpp )
target_link_libraries ( MESA_PD_Kernel_CNT_WallContact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CNT_WallContact )

waLBerla_add_test_executable ( MESA_PD_Kernel_CoefficientOfRestitutionSD kernel/CoefficientOfRestitutionSD.cpp )
target_link_libraries ( MESA_PD_Kernel_CoefficientOfRestitutionSD walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionSDVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionSD> --useVV )

waLBerla_add_test_executable ( MESA_PD_Kernel_CoefficientOfRestitutionLSD kernel/CoefficientOfRestitutionLSD.cpp )
target_link_libraries ( MESA_PD_Kernel_CoefficientOfRestitutionLSD walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionLSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionLSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionLSDVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionLSD> --useVV )

waLBerla_add_test_executable ( MESA_PD_Kernel_CoefficientOfRestitutionNLSD kernel/CoefficientOfRestitutionNLSD.cpp )
target_link_libraries ( MESA_PD_Kernel_CoefficientOfRestitutionNLSD walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionNLSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionNLSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionNLSD> --useVV )

waLBerla_add_test_executable ( MESA_PD_Kernel_DetectAndStoreContacts kernel/DetectAndStoreContacts.cpp )
target_link_libraries ( MESA_PD_Kernel_DetectAndStoreContacts walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_DetectAndStoreContacts )

waLBerla_add_test_executable ( MESA_PD_Kernel_DoubleCast kernel/DoubleCast.cpp )
target_link_libraries ( MESA_PD_Kernel_DoubleCast walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_DoubleCast )

waLBerla_add_test_executable ( MESA_PD_Kernel_ExplicitEuler kernel/ExplicitEuler.cpp )
target_link_libraries ( MESA_PD_Kernel_ExplicitEuler walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ExplicitEuler )

waLBerla_add_test_executable ( MESA_PD_Kernel_ForceLJ kernel/ForceLJ.cpp )
target_link_libraries ( MESA_PD_Kernel_ForceLJ walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ForceLJ )

waLBerla_add_test_executable ( MESA_PD_Kernel_GenerateAnalyticContacts kernel/GenerateAnalyticContacts.cpp )
target_link_libraries ( MESA_PD_Kernel_GenerateAnalyticContacts walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_GenerateAnalyticContacts PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_GenerateLinkedCells kernel/GenerateLinkedCells.cpp )
target_link_libraries ( MESA_PD_Kernel_GenerateLinkedCells walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_GenerateLinkedCells )

waLBerla_add_test_executable ( MESA_PD_Kernel_HCSITSKernels kernel/HCSITSKernels.cpp )
target_link_libraries ( MESA_PD_Kernel_HCSITSKernels walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_HCSITSKernels )

waLBerla_add_test_executable ( MESA_PD_Kernel_HeatConduction kernel/HeatConduction.cpp )
target_link_libraries ( MESA_PD_Kernel_HeatConduction walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_HeatConduction )

waLBerla_add_test_executable ( MESA_PD_Kernel_IntegratorAccuracy kernel/IntegratorAccuracy.cpp )
target_link_libraries ( MESA_PD_Kernel_IntegratorAccuracy walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_IntegratorAccuracyEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_IntegratorAccuracy> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_IntegratorAccuracyVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_IntegratorAccuracy> --useVV )

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/kernel/interfaces MESA_PD_INTERFACE_CHECKS)
waLBerla_add_test_executable ( MESA_PD_Kernel_Interfaces kernel/Interfaces.cpp ${MESA_PD_INTERFACE_CHECKS} )
target_link_libraries ( MESA_PD_Kernel_Interfaces walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_Interfaces )

waLBerla_add_test_executable ( MESA_PD_Kernel_LinearSpringDashpot kernel/LinearSpringDashpot.cpp )
target_link_libraries ( MESA_PD_Kernel_LinearSpringDashpot walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_LinearSpringDashpot )

waLBerla_add_test_executable ( MESA_PD_Kernel_LinkedCellsVsBruteForce kernel/LinkedCellsVsBruteForce.cpp )
target_link_libraries ( MESA_PD_Kernel_LinkedCellsVsBruteForce walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_LinkedCellsVsBruteForce PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_PFCDamping kernel/PFCDamping.cpp )
target_link_libraries ( MESA_PD_Kernel_PFCDamping walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_PFCDamping )

waLBerla_add_test_executable ( MESA_PD_Kernel_SemiImplicitEuler kernel/SemiImplicitEuler.cpp )
target_link_libraries ( MESA_PD_Kernel_SemiImplicitEuler walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SemiImplicitEuler )

waLBerla_add_test_executable ( MESA_PD_Kernel_SingleCast kernel/SingleCast.cpp )
target_link_libraries ( MESA_PD_Kernel_SingleCast walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SingleCast )

waLBerla_add_test_executable ( MESA_PD_Kernel_SpherePile kernel/SpherePile.cpp )
target_link_libraries ( MESA_PD_Kernel_SpherePile walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SpherePile )

waLBerla_add_test_executable ( MESA_PD_Kernel_SpringDashpot kernel/SpringDashpot.cpp )
target_link_libraries ( MESA_PD_Kernel_SpringDashpot walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SpringDashpot )

waLBerla_add_test_executable ( MESA_PD_Kernel_SyncGhostOwners kernel/SyncGhostOwners.cpp )
target_link_libraries ( MESA_PD_Kernel_SyncGhostOwners walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncGhostOwners PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_SyncGhostOwnersLarge kernel/SyncGhostOwnersLarge.cpp )
target_link_libraries ( MESA_PD_Kernel_SyncGhostOwnersLarge walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncGhostOwnersLarge PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_SyncNextNeighbors kernel/SyncNextNeighbors.cpp )
target_link_libraries ( MESA_PD_Kernel_SyncNextNeighbors walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncNextNeighbors PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_SyncNextNeighborsBlockForest kernel/SyncNextNeighborsBlockForest.cpp )
target_link_libraries ( MESA_PD_Kernel_SyncNextNeighborsBlockForest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncNextNeighborsBlockForest PROCESSES 27 )

waLBerla_add_test_executable ( MESA_PD_Kernel_TemperatureIntegration kernel/TemperatureIntegration.cpp )
target_link_libraries ( MESA_PD_Kernel_TemperatureIntegration walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_TemperatureIntegration )

waLBerla_add_test_executable ( MESA_PD_Kernel_VelocityVerlet kernel/VelocityVerlet.cpp )
target_link_libraries ( MESA_PD_Kernel_VelocityVerlet walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Kernel_VelocityVerlet )

waLBerla_add_test_executable ( MESA_PD_MPI_BroadcastProperty mpi/BroadcastProperty.cpp )
target_link_libraries ( MESA_PD_MPI_BroadcastProperty walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_BroadcastProperty PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_MPI_ClearGhostOwnerSync mpi/ClearGhostOwnerSync.cpp )
target_link_libraries ( MESA_PD_MPI_ClearGhostOwnerSync walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_ClearGhostOwnerSync PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_MPI_ClearNextNeighborSync mpi/ClearNextNeighborSync.cpp )
target_link_libraries ( MESA_PD_MPI_ClearNextNeighborSync walberla::blockforest walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_ClearNextNeighborSync PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_MPI_Notifications mpi/Notifications.cpp )
target_link_libraries ( MESA_PD_MPI_Notifications walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_Notifications )

waLBerla_add_test_executable ( MESA_PD_MPI_ReduceContactHistory mpi/ReduceContactHistory.cpp )
target_link_libraries ( MESA_PD_MPI_ReduceContactHistory walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_ReduceContactHistory PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_MPI_ReduceProperty mpi/ReduceProperty.cpp )
target_link_libraries ( MESA_PD_MPI_ReduceProperty walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_ReduceProperty PROCESSES 8 )

waLBerla_add_test_executable ( MESA_PD_MPI_ShapePackUnpack mpi/ShapePackUnpack.cpp )
target_link_libraries ( MESA_PD_MPI_ShapePackUnpack walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_ShapePackUnpack )

waLBerla_add_test_executable ( MESA_PD_MPI_SyncNextNeighborsNoGhosts mpi/SyncNextNeighborsNoGhosts.cpp )
target_link_libraries ( MESA_PD_MPI_SyncNextNeighborsNoGhosts walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_SyncNextNeighborsNoGhosts PROCESSES 2 )

waLBerla_add_test_executable ( MESA_PD_MPI_VelocityCorrectionNotification mpi/VelocityCorrectionNotification.cpp )
target_link_libraries ( MESA_PD_MPI_VelocityCorrectionNotification walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_MPI_VelocityCorrectionNotification PROCESSES 8)

waLBerla_add_test_executable ( MESA_PD_Sorting Sorting.cpp )
target_link_libraries ( MESA_PD_Sorting walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_Sorting )

waLBerla_add_test_executable ( MESA_PD_Stiffness Stiffness.cpp )
target_link_libraries ( MESA_PD_Stiffness walberla::blockforest walberla::core walberla::mesa_pd walberla::mesa_pd )

waLBerla_add_test_executable ( MESA_PD_VTK_Outputs vtk/VTKOutputs.cpp )
target_link_libraries ( MESA_PD_VTK_Outputs walberla::blockforest walberla::core walberla::vtk walberla::mesa_pd )
waLBerla_execute_test( NAME   MESA_PD_VTK_Outputs PROCESSES 8 )

if (WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
    waLBerla_add_test_executable ( MESA_PD_VTK_ConvexPolyhedron vtk/ConvexPolyhedronVTKOutput.cpp )
    target_link_libraries ( MESA_PD_VTK_ConvexPolyhedron walberla::core walberla::mesa_pd walberla::mesh_common walberla::vtk walberla::mesa_pd )
    waLBerla_execute_test( NAME   MESA_PD_VTK_ConvexPolyhedron )

    waLBerla_add_test_executable ( MESA_PD_Common_ContainsPointConvexPolyhedron common/ContainsPointConvexPolyhedron.cpp )
    target_link_libraries ( MESA_PD_Common_ContainsPointConvexPolyhedron walberla::core walberla::mesa_pd walberla::mesh_common walberla::mesa_pd )
    waLBerla_execute_test( NAME   MESA_PD_Common_ContainsPointConvexPolyhedron )

    waLBerla_add_test_executable ( MESA_PD_Data_ConvexPolyhedron data/ConvexPolyhedron.cpp )
    target_link_libraries ( MESA_PD_Data_ConvexPolyhedron walberla::core walberla::mesa_pd walberla::mesh_common walberla::mesa_pd )
    waLBerla_execute_test( NAME   MESA_PD_Data_ConvexPolyhedron )

    waLBerla_add_test_executable ( MESA_PD_COLLISIONDETECTION_ConvexPolyhedron_GJK_EPA collision_detection/ConvexPolyhedron_GJK_EPA.cpp )
    target_link_libraries ( MESA_PD_COLLISIONDETECTION_ConvexPolyhedron_GJK_EPA walberla::core walberla::mesa_pd walberla::mesh_common walberla::mesa_pd )
    waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_ConvexPolyhedron_GJK_EPA )
endif()
