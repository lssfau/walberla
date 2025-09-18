###################################################################################################
#
# Module executiontree
#
###################################################################################################

add_library( walberla_executiontree )
add_library( walberla::executiontree ALIAS walberla_executiontree )
target_link_libraries( walberla_executiontree PUBLIC walberla::core walberla::domain_decomposition walberla::timeloop )

target_sources( walberla_executiontree
        PUBLIC
        ExecutionTreeSweep.h
        ExecutionTree.h
        ExecutionTree.impl.h
)

target_sources( walberla_executiontree
        PRIVATE
        ExecutionTree.cpp
)
                        
###################################################################################################                        