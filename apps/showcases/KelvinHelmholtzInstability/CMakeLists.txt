waLBerla_link_files_to_builddir( *.py )

if( WALBERLA_BUILD_WITH_CODEGEN )
    foreach(collision_setup cumulant-allOne cumulant cumulant-L1p000 cumulant-L0p100 cumulant-L0p010 cumulant-L0p001 smagorinsky qr)
        set(config ${collision_setup})

        walberla_generate_target_from_python( 
            NAME KHI_CodeGeneration_${config}
            FILE KHI.py
            CODEGEN_CFG ${config}
            OUT_FILES 
                KHIStorageSpecification.h KHIStorageSpecification.${CODEGEN_FILE_SUFFIX}
                KHISweepCollection.h KHISweepCollection.${CODEGEN_FILE_SUFFIX}
                KHIBoundaryCollection.h
                KHI_AdvectionDiffusionStorageSpecification.h KHI_AdvectionDiffusionStorageSpecification.${CODEGEN_FILE_SUFFIX}
                KHI_AdvectionDiffusionSweepCollection.h KHI_AdvectionDiffusionSweepCollection.${CODEGEN_FILE_SUFFIX}
                KHI_AdvectionDiffusionBoundaryCollection.h
                FreeSlip_BC.h FreeSlip_BC.${CODEGEN_FILE_SUFFIX}
                Neumann_BC.h Neumann_BC.${CODEGEN_FILE_SUFFIX}
                InfoHeader.h
        )

        set(KHI_DEPENDENCIES walberla::blockforest walberla::core walberla::domain_decomposition walberla::field walberla::geometry walberla::timeloop walberla::lbm walberla::python_coupling walberla::stencil walberla::vtk KHI_CodeGeneration_${config})
    
        if(WALBERLA_BUILD_WITH_GPU_SUPPORT)
            list(APPEND KHI_DEPENDENCIES walberla::gpu)
        endif()

        walberla_add_executable( 
            NAME KHI_${config}
            FILES KHI.cpp
            DEPENDS ${KHI_DEPENDENCIES}
        )
    endforeach()
endif()
