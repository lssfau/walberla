message(STATUS "Setting MSVC specific compiler options")

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.11)
  message(FATAL_ERROR "MSVC version must be at least 19.11!")
endif()

if(WALBERLA_PROFILE_GENERATE)
  add_flag(CMAKE_CXX_FLAGS "/GL")
  add_flag(CMAKE_MODULE_LINKER_FLAGS "/LTCG:PGINSTRUMENT")
  add_flag(CMAKE_SHARED_LINKER_FLAGS "/LTCG:PGINSTRUMENT")
  add_flag(CMAKE_EXE_LINKER_FLAGS "/LTCG:PGINSTRUMENT")
endif()

if(WALBERLA_PROFILE_USE)
  add_flag(CMAKE_CXX_FLAGS "/GL")
  add_flag(CMAKE_MODULE_LINKER_FLAGS "/LTCG:PGOPTIMIZE")
  add_flag(CMAKE_SHARED_LINKER_FLAGS "/LTCG:PGOPTIMIZE")
  add_flag(CMAKE_EXE_LINKER_FLAGS "/LTCG:PGOPTIMIZE")
endif()

if(WALBERLA_BUILD_WITH_FASTMATH)
  add_flag(CMAKE_CXX_FLAGS "/fp:fast")
endif()

string(REGEX REPLACE "[/-]W[0-4]" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
)# remove default warning flags

option(WALBERLA_GROUP_PROJECTS
       "Flag if the projects are grouped or in a flat hierarchy" ON)
option(WALBERLA_GROUP_FILES
       "Flag if the files are grouped or in a flat hierarchy" ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ${WALBERLA_GROUP_PROJECTS})

option(WALBERLA_VS_MULTI_PROCESS_BUILD "Use the /mp option for VS builds" ON)
if(WALBERLA_VS_MULTI_PROCESS_BUILD)
  add_flag(CMAKE_CXX_FLAGS "-MP") # enable multi-threaded compiling
endif()

add_definitions("-DNOMINMAX") # Disable Min/Max-Macros
add_definitions("-D_WIN32_WINNT=0x501") # Minimum Windows versions is Windows XP
add_definitions("-DWINVER=0x501") # Minimum Windows versions is Windows XP
add_definitions("-D_CRT_SECURE_NO_WARNINGS") # disable warnings promoting
                                             # Microsoft's security enhanced CRT
add_definitions("-D_SCL_SECURE_NO_WARNINGS") # disable warnings triggered by
                                             # Microsoft's checked iterators
add_flag(CMAKE_CXX_FLAGS "-W4") # set warning level to maximum
add_flag(CMAKE_CXX_FLAGS "-bigobj") # enable big object files
add_flag(CMAKE_CXX_FLAGS "-wd4127") # disable compiler warning C4127:
                                    # "conditional expression is constant"
add_flag(CMAKE_CXX_FLAGS "-wd4512") # disable compiler warning C4512:
                                    # "assignment operator could not be
                                    # generated"
add_flag(CMAKE_CXX_FLAGS "-wd4913") # disable compiler warning C4512: "user
                                    # defined binary operator ',' exists but
# no overload could convert all operands, default built-in binary operator ','
# used"
add_flag(CMAKE_CXX_FLAGS "-wd4702") # disable compiler warning C4702:
                                    # "unreachable code"
add_flag(CMAKE_CXX_FLAGS "-wd4505") # disable compiler warning C4505:
                                    # "unreferenced local function has been
                                    # removed"
add_flag(CMAKE_CXX_FLAGS "-wd4503") # disable compiler warning C4503:
                                    # "'identifier' : decorated name length
                                    # exceeded, name was truncated"

if(NOT WARNING_DEPRECATED)
  add_definitions("-D_CRT_SECURE_NO_DEPRECATE")
  add_definitions("-D_SCL_SECURE_NO_DEPRECATE")
  add_flag(CMAKE_CXX_FLAGS "-wd4996") # Disable compiler warning C4996:
                                      # "declared as deprecated"
endif()

string(REPLACE "/Od" "/O2" CMAKE_C_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
string(REPLACE "/Ob0" "/Ob2" CMAKE_C_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
string(REPLACE "/RTC1" "" CMAKE_C_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
string(REPLACE "/Od" "/O2" CMAKE_CXX_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
string(REPLACE "/Ob0" "/Ob2" CMAKE_CXX_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
string(REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUGOPTIMIZED
               ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
