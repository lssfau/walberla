add_subdirectory(setups)

waLBerla_add_executable(NAME SettlingSpheres
        FILES SettlingSpheres.cpp
      DEPENDS walberla::blockforest walberla::core walberla::field walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::vtk )
