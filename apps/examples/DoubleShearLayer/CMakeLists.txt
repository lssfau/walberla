walberla_link_files_to_builddir( *.prm )

add_executable( DoubleShearLayer )
target_sources( DoubleShearLayer PRIVATE DoubleShearLayer.cpp )

walberla_generate_sources(
    DoubleShearLayer
    SCRIPTS DoubleShearLayerSweeps.py
)

target_link_libraries(
    DoubleShearLayer
    PRIVATE
    walberla::core
    walberla::experimental
    walberla::blockforest
    walberla::geometry
    walberla::stencil
    walberla::vtk
)
