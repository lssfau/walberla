walberla_link_files_to_builddir( *.prm )

add_executable( DoubleShearLayer )
target_sources( DoubleShearLayer PRIVATE DoubleShearLayer.cpp )
walberla_set_gpu_language( DoubleShearLayer.cpp )

walberla_generate_sources(
    DoubleShearLayer
    SCRIPTS DoubleShearLayerSweeps.py
    AUTO_LANGUAGE
)

target_link_libraries( DoubleShearLayer PRIVATE walberla::v8 )
