###################################################################################################
#
# Module postprocessing
#
###################################################################################################

add_library( walberla_postprocessing )
add_library( walberla::postprocessing ALIAS walberla_postprocessing )
target_link_libraries( walberla_postprocessing PUBLIC walberla::core walberla::domain_decomposition walberla::field walberla::geometry )
target_sources( walberla_postprocessing
        PUBLIC
        MarchingCubes.impl.h
        MarchingCubes.h
        FieldToSurfaceMesh.h
)

target_sources( walberla_postprocessing
        PRIVATE
        FieldToSurfaceMesh.cpp
)

###################################################################################################
