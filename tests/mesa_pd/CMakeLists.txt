###################################################################################################
#
# Tests for MESA_PD
#
###################################################################################################

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticCollisionFunctions FILES collision_detection/AnalyticCollisionFunctions.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticCollisionFunctions )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticContactDetection FILES collision_detection/AnalyticContactDetection.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_AnalyticContactDetection )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_BoxSupport FILES collision_detection/BoxSupport.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_BoxSupport )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_EllipsoidSupport FILES collision_detection/EllipsoidSupport.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_EllipsoidSupport )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_EPA FILES collision_detection/EPA.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_EPA )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_GeneralContactDetection FILES collision_detection/GeneralContactDetection.cpp DEPENDS blockforest core pe)
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GeneralContactDetection )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_GJK FILES collision_detection/GJK.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GJK )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_GJKEPA FILES collision_detection/GJK_EPA.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_GJKEPA )

waLBerla_compile_test( NAME   MESA_PD_COLLISIONDETECTION_SphereSupport FILES collision_detection/SphereSupport.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COLLISIONDETECTION_SphereSupport )

waLBerla_compile_test( NAME   MESA_PD_COMMON_IntersectionRatio FILES common/IntersectionRatio.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_COMMON_IntersectionRatio )

waLBerla_compile_test( NAME   MESA_PD_ContactDetection FILES ContactDetection.cpp DEPENDS blockforest core pe)
waLBerla_execute_test( NAME   MESA_PD_ContactDetection PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_Data_Flags FILES data/Flags.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Data_Flags )

waLBerla_compile_test( NAME   MESA_PD_Data_ParticleStorage FILES data/ParticleStorage.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Data_ParticleStorage )

waLBerla_compile_test( NAME   MESA_PD_Domain_BlockForestDomain FILES domain/BlockForestDomain.cpp DEPENDS blockforest core )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestDomain )

waLBerla_compile_test( NAME   MESA_PD_Domain_BlockForestSync FILES domain/BlockForestSync.cpp DEPENDS blockforest core pe)
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync2 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 2 )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync4 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 4 )
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSync8 COMMAND $<TARGET_FILE:MESA_PD_Domain_BlockForestSync> PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_Domain_BlockForestSyncPeriodic FILES domain/BlockForestSyncPeriodic.cpp DEPENDS blockforest core pe)
waLBerla_execute_test( NAME   MESA_PD_Domain_BlockForestSyncPeriodic PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_Domain_DistanceCalculation FILES domain/DistanceCalculation.cpp DEPENDS blockforest core )
waLBerla_execute_test( NAME   MESA_PD_Domain_DistanceCalculation )

waLBerla_compile_test( NAME   MESA_PD_Domain_DynamicRefinement FILES domain/DynamicRefinement.cpp DEPENDS blockforest core pe )
waLBerla_execute_test( NAME   MESA_PD_Domain_DynamicRefinement PROCESSES 8)

waLBerla_compile_test( NAME   MESA_PD_Domain_SerializeDeserialize FILES domain/SerializeDeserialize.cpp DEPENDS blockforest core pe)
waLBerla_execute_test( NAME   MESA_PD_Domain_SerializeDeserialize PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_DropTestAnalytic FILES DropTestAnalytic.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_DropTestAnalytic )

waLBerla_compile_test( NAME   MESA_PD_DropTestGeneral FILES DropTestGeneral.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_DropTestGeneral )

waLBerla_compile_test( NAME   MESA_PD_Kernel_ClearNextNeighborSync FILES kernel/ClearNextNeighborSync.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ClearNextNeighborSync PROCESSES 2 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionSD FILES kernel/CoefficientOfRestitutionSD.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionSDVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionSD> --useVV )

waLBerla_compile_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionLSD FILES kernel/CoefficientOfRestitutionLSD.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionLSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionLSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionLSDVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionLSD> --useVV )

waLBerla_compile_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionNLSD FILES kernel/CoefficientOfRestitutionNLSD.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionNLSDEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionNLSD> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_CoefficientOfRestitutionVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_CoefficientOfRestitutionNLSD> --useVV )

waLBerla_compile_test( NAME   MESA_PD_Kernel_DetectAndStoreContacts FILES kernel/DetectAndStoreContacts.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_DetectAndStoreContacts )

waLBerla_compile_test( NAME   MESA_PD_Kernel_DoubleCast FILES kernel/DoubleCast.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_DoubleCast )

waLBerla_compile_test( NAME   MESA_PD_Kernel_ExplicitEuler FILES kernel/ExplicitEuler.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ExplicitEuler )

waLBerla_compile_test( NAME   MESA_PD_Kernel_ExplicitEulerWithShape FILES kernel/ExplicitEulerWithShape.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ExplicitEulerWithShape )

waLBerla_compile_test( NAME   MESA_PD_Kernel_ForceLJ FILES kernel/ForceLJ.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_ForceLJ )

waLBerla_compile_test( NAME   MESA_PD_Kernel_GenerateAnalyticContacts FILES kernel/GenerateAnalyticContacts.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_GenerateAnalyticContacts PROCESSES 27 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_GenerateLinkedCells FILES kernel/GenerateLinkedCells.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_GenerateLinkedCells )

waLBerla_compile_test( NAME   MESA_PD_Kernel_HeatConduction FILES kernel/HeatConduction.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_HeatConduction )

waLBerla_compile_test( NAME   MESA_PD_Kernel_IntegratorAccuracy FILES kernel/IntegratorAccuracy.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_IntegratorAccuracyEuler COMMAND $<TARGET_FILE:MESA_PD_Kernel_IntegratorAccuracy> )
waLBerla_execute_test( NAME   MESA_PD_Kernel_IntegratorAccuracyVelocityVerlet COMMAND $<TARGET_FILE:MESA_PD_Kernel_IntegratorAccuracy> --useVV )

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/kernel/interfaces MESA_PD_INTERFACE_CHECKS)
waLBerla_compile_test( NAME   MESA_PD_Kernel_Interfaces FILES kernel/Interfaces.cpp ${MESA_PD_INTERFACE_CHECKS} DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_Interfaces )

waLBerla_compile_test( NAME   MESA_PD_Kernel_LinearSpringDashpot FILES kernel/LinearSpringDashpot.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_LinearSpringDashpot )

waLBerla_compile_test( NAME   MESA_PD_Kernel_LinkedCellsVsBruteForce FILES kernel/LinkedCellsVsBruteForce.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_LinkedCellsVsBruteForce PROCESSES 27 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_SingleCast FILES kernel/SingleCast.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SingleCast )

waLBerla_compile_test( NAME   MESA_PD_Kernel_SpringDashpot FILES kernel/SpringDashpot.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SpringDashpot )

waLBerla_compile_test( NAME   MESA_PD_Kernel_SyncGhostOwners FILES kernel/SyncGhostOwners.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncGhostOwners PROCESSES 27 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_SyncGhostOwnersLarge FILES kernel/SyncGhostOwnersLarge.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncGhostOwnersLarge PROCESSES 27 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_SyncNextNeighbors FILES kernel/SyncNextNeighbors.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_SyncNextNeighbors PROCESSES 27 )

waLBerla_compile_test( NAME   MESA_PD_Kernel_TemperatureIntegration FILES kernel/TemperatureIntegration.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_TemperatureIntegration )

waLBerla_compile_test( NAME   MESA_PD_Kernel_VelocityVerlet FILES kernel/VelocityVerlet.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_VelocityVerlet )

waLBerla_compile_test( NAME   MESA_PD_Kernel_VelocityVerletWithShape FILES kernel/VelocityVerletWithShape.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Kernel_VelocityVerletWithShape )

waLBerla_compile_test( NAME   MESA_PD_MPI_BroadcastProperty FILES mpi/BroadcastProperty.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_MPI_BroadcastProperty PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_MPI_Notifications FILES mpi/Notifications.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_MPI_Notifications )

waLBerla_compile_test( NAME   MESA_PD_MPI_ReduceContactHistory FILES mpi/ReduceContactHistory.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_MPI_ReduceContactHistory PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_MPI_ReduceProperty FILES mpi/ReduceProperty.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_MPI_ReduceProperty PROCESSES 8 )

waLBerla_compile_test( NAME   MESA_PD_Sorting FILES Sorting.cpp DEPENDS core )
waLBerla_execute_test( NAME   MESA_PD_Sorting )

waLBerla_compile_test( NAME   MESA_PD_VTK_Outputs FILES vtk/VTKOutputs.cpp DEPENDS blockforest core vtk )
waLBerla_execute_test( NAME   MESA_PD_VTK_Outputs PROCESSES 8 )
