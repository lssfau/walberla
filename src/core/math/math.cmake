
OPTION( WALBERLA_DEACTIVATE_MATH_PARSER "Deactivates the function parser of core/math" OFF )

mark_as_advanced( WALBERLA_DEACTIVATE_MATH_PARSER )

configure_file ( math/CMakeDefs.in.h math/CMakeDefs.h )
