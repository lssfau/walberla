mkdir build
cd build

cmake \
   -DCMAKE_FIND_ROOT_PATH=${PREFIX} \
   -DCMAKE_INSTALL_PREFIX=${PREFIX} \
   -DWALBERLA_BUILD_WITH_PYTHON=ON \
   ..

make -j ${CPU_COUNT} pythonModuleInstall
