waLBerla_link_files_to_builddir( *.prm )

if( WALBERLA_BUILD_WITH_CODEGEN )

    foreach(stencil d3q19 d3q27)
        foreach(weights none lattice geometric)
            foreach(velocity filtered mean)
                foreach(sgs_model none smago qr)

                    if (${stencil} STREQUAL "d3q19" AND ( NOT ${weights} STREQUAL "none" ))
                        # D3Q19 stencil does not need a closure for the weights as it is complete
                        continue()
                    endif()

                    if (${stencil} STREQUAL "d3q27" AND ${weights} STREQUAL "none" )
                        # D3Q27 stencil needs a closure for the weights as it is incomplete
                        continue()
                    endif()

                    set(config ${stencil}_${weights}_${velocity}_${sgs_model})

                    #   Turbulent Channel generation
                    walberla_generate_target_from_python( NAME TurbulentChannel_CodeGeneration_${config}
                            FILE TurbulentChannel.py
                            CODEGEN_CFG ${config}
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

                    add_executable( TurbulentChannel_Application_${config} TurbulentChannel.cpp )
                    target_link_libraries( TurbulentChannel_Application_${config}
                            walberla::blockforest walberla::core walberla::domain_decomposition walberla::field
                            walberla::geometry walberla::timeloop walberla::lbm walberla::stencil walberla::vtk
                            TurbulentChannel_CodeGeneration_${config} )
                endforeach ()
            endforeach ()
        endforeach ()
    endforeach()

endif()