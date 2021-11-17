
mkdir build
cd build

set CMAKE_CONFIG="Release"
set PYTHON_LIBRARY=%PREFIX%\libs\python%PY_VER:~0,1%%PY_VER:~2,1%.lib

cmake -LAH -G"Visual Studio 15 2017 Win64"                   ^
  -DCMAKE_BUILD_TYPE="%CMAKE_CONFIG%"                        ^
  -DCMAKE_FIND_ROOT_PATH="%LIBRARY_PREFIX%"                  ^
  -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%"                  ^
  -DBUILD_APPS=OFF                                           ^
  -DPython_EXECUTABLE:FILEPATH="%PYTHON_LIBRARY%"               ^
  -DOPENMESH_PYTHON_VERSION="%PY_VER%"                       ^
  -DPYTHON_INSTALL_DIR="%SP_DIR%"                            ^
  -DOPENMESH_BUILD_PYTHON_UNIT_TESTS=ON ..
if errorlevel 1 exit 1

cmake --build . --config %CMAKE_CONFIG% --target INSTALL
if errorlevel 1 exit 1

move %LIBRARY_PREFIX%\lib\python\openmesh.lib %SP_DIR%
move %LIBRARY_PREFIX%\lib\python\openmesh.pyd %SP_DIR%

ctest -C %CMAKE_CONFIG% --output-on-failure
if errorlevel 1 exit 1
