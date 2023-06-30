###################################################################################################
#
# Module timeloop
#
###################################################################################################

add_library( timeloop )
target_link_libraries( timeloop PUBLIC core domain_decomposition )
target_sources( timeloop
      PRIVATE
      ITimeloop.h
      PerformanceMeter.h
      SelectableFunctionCreators.h
      SweepTimeloop.h
      Timeloop.h
      PerformanceMeter.cpp
)

###################################################################################################
