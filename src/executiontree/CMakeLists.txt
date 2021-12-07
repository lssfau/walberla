###################################################################################################
#
# Module executiontree
#
###################################################################################################

add_library( executiontree )
target_link_libraries( executiontree PUBLIC core domain_decomposition timeloop )
target_sources( executiontree
      PRIVATE
      ExecutionTreeSweep.h
      ExecutionTree.cpp
      ExecutionTree.h
      ExecutionTree.impl.h
      )
                        
###################################################################################################                        