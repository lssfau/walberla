add_subdirectory(setups)

waLBerla_add_executable(NAME SettlingSpheres
        FILES SettlingSpheres.cpp
        DEPENDS blockforest core field lbm_mesapd_coupling mesa_pd vtk)
