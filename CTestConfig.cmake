## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "waLBerla")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "www10.informatik.uni-erlangen.de")
set(CTEST_DROP_LOCATION "/Research/Projects/walberla/CDash/submit.php?project=waLBerla")
set(CTEST_DROP_SITE_CDASH TRUE)

set (MEMORYCHECK_COMMAND:FILEPATH "/usr/bin/valgrind")
set (MEMORYCHECK_COMMAND_OPTIONS "--show-reachable=no")
set (MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/openmpi-valgrind.supp")
