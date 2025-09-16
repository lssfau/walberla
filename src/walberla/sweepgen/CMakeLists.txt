add_library( walberla_sweepgen INTERFACE )
add_library( walberla::sweepgen ALIAS walberla_sweepgen )

target_sources( walberla_experimental
    INTERFACE
    boundaries/GenericHbbBoundary.hpp
    boundaries/IrregularFreeSlip.hpp
    communication/UniformGpuFieldPackInfoBase.hpp
)

target_link_libraries(
    walberla_sweepgen
    INTERFACE
    walberla_core
    walberla_stencil
    walberla_blockforest
    walberla_field
)

if( WALBERLA_BUILD_WITH_GPU )
    target_link_libraries( walberla_sweepgen INTERFACE walberla_gpu )
endif()