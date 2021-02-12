###################################################################################################
#
# Tests for timeloop module
#
###################################################################################################


waLBerla_compile_test( FILES TimeloopAndSweepRegister.cpp DEPENDS field blockforest )
waLBerla_execute_test(NAME TimeloopAndSweepRegister )
