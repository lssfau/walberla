###################################################################################################
#
# Module timeloop
#
###################################################################################################

add_library( walberla_timeloop )
add_library( walberla::timeloop ALIAS walberla_timeloop )
target_link_libraries( walberla_timeloop PUBLIC walberla::core walberla::domain_decomposition )
target_sources( walberla_timeloop
      PRIVATE
      ITimeloop.h
      PerformanceMeter.h
      SelectableFunctionCreators.h
      SweepTimeloop.h
      Timeloop.h
      PerformanceMeter.cpp
)

###################################################################################################
