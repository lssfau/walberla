
OPTION( WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID "BlockForest: enables/disables BlockID implementation as a primitive data type (uint_t)" ON )

mark_as_advanced( WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID )

configure_file ( CMakeDefs.in.h  CMakeDefs.h )

waLBerla_add_module( DEPENDS communication core domain_decomposition stencil )
