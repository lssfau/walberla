waLBerla_add_executable ( NAME VBondModel
                          FILES VBondModel.cpp
                          DEPENDS core mesa_pd )
#waLBerla_execute_test( NAME MESA_PD_CNT_VBondModel
#                       COMMAND $<TARGET_FILE:VBondModel> )
