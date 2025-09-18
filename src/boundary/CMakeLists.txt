###################################################################################################
#
# Module boundary
#
###################################################################################################


add_library( walberla_boundary )
add_library( walberla::boundary ALIAS walberla_boundary )
target_link_libraries( walberla_boundary PUBLIC walberla::core walberla::communication walberla::domain_decomposition walberla::field walberla::stencil )
target_sources( walberla_boundary
        PUBLIC
        BoundaryHandling.h
        Boundary.h
        BoundaryHandlingCollection.h
        BoundaryUID.h
        ShiftedPeriodicity.h
)

target_sources( walberla_boundary
        PRIVATE
        Boundary.cpp
)
