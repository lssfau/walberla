mkdir build
cd build

set BOOST_ROOT=%PREFIX%
cmake -LAH -G"Visual Studio 15 2017 Win64"                   ^
  -DWALBERLA_BUILD_WITH_PYTHON=ON                            ^
  -DWALBERLA_BUILD_WITH_MPI=OFF                              ^
  -DWALBERLA_BUILD_WITH_OPENMP=ON                            ^
  -DPYTHON_EXECUTABLE="%PYTHON%"  ..
if errorlevel 1 exit 1

cmake --build . --config Release --target pythonModuleInstall
if errorlevel 1 exit 1