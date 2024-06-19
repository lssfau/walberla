waLBerla_link_files_to_builddir( *.prm )

if( WALBERLA_BUILD_WITH_CODEGEN )

    #   Turbulent Channel generation
    walberla_generate_target_from_python( NAME TurbulentChannel_CodeGeneration
            FILE TurbulentChannel.py
            OUT_FILES CodegenIncludes.h
            TurbulentChannel_Sweep.cpp TurbulentChannel_Sweep.h
            TurbulentChannel_PackInfo.cpp TurbulentChannel_PackInfo.h
            TurbulentChannel_Setter.cpp TurbulentChannel_Setter.h
            TurbulentChannel_NoSlip.cpp TurbulentChannel_NoSlip.h
            TurbulentChannel_FreeSlip_top.cpp TurbulentChannel_FreeSlip_top.h
            TurbulentChannel_WFB_bottom.cpp TurbulentChannel_WFB_bottom.h
            TurbulentChannel_WFB_top.cpp TurbulentChannel_WFB_top.h
            TurbulentChannel_Welford.cpp TurbulentChannel_Welford.h
            TurbulentChannel_Welford_TKE_SGS.cpp TurbulentChannel_Welford_TKE_SGS.h
            TurbulentChannel_TKE_SGS_Writer.cpp TurbulentChannel_TKE_SGS_Writer.h
    )

    walberla_add_executable ( NAME TurbulentChannel_Application
            FILES TurbulentChannel.cpp
            DEPENDS blockforest core domain_decomposition field geometry timeloop lbm stencil vtk
                    TurbulentChannel_CodeGeneration )

endif()