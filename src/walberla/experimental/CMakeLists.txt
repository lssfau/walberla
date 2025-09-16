add_library( walberla_experimental INTERFACE )
add_library( walberla::experimental ALIAS walberla_experimental )

target_sources( walberla_experimental
    INTERFACE
    memory/MemoryTags.hpp
    memory/UnifiedMemoryAllocator.hpp
    sweep/SparseIndexList.hpp
)

target_link_libraries(
    walberla_experimental
    INTERFACE
    walberla_core
    walberla_stencil
    walberla_blockforest
    walberla_field
)

if( WALBERLA_BUILD_WITH_GPU )
    target_link_libraries( walberla_experimental INTERFACE walberla_gpu )
endif()
