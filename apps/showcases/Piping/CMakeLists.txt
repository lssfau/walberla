add_subdirectory(setups)

add_executable ( SettlingSpheres SettlingSpheres.cpp )
target_link_libraries ( SettlingSpheres walberla::blockforest walberla::core walberla::field
        walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::vtk )
