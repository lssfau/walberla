waLBerla_link_files_to_builddir( *.dat )

add_executable ( 01_cnt_film 01_cnt_film.cpp Config.cpp FilmSpecimen.cpp InitializeCNTs.cpp SQLProperties.cpp
        Statistics.cpp )

target_link_libraries ( 01_cnt_film walberla::blockforest walberla::core walberla::mesa_pd walberla::sqlite walberla::vtk )
