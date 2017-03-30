from distutils.core import setup
import shutil
from os.path import exists, join

# The following variables are configure by CMake
walberla_source_dir  = "${walberla_SOURCE_DIR}"
walberla_binary_dir = "${CMAKE_CURRENT_BINARY_DIR}"

def collectFiles():
    shared_lib = join(walberla_binary_dir, 'walberla_cpp.so' )
    # copy everything inplace
    if not exists( shared_lib ):
        print("Python Module was not built yet - run 'make walberla_cpp'")
        exit(1)


    shutil.rmtree( join(walberla_binary_dir, 'waLBerla'), ignore_errors=True )

    shutil.copytree( join(walberla_source_dir, 'python', 'waLBerla'),
                     join(walberla_binary_dir, 'waLBerla') )

    shutil.copy( shared_lib,
                 join(walberla_binary_dir, 'waLBerla', 'walberla_cpp.so') )

packages = ['waLBerla',
            'waLBerla.evaluation',
            'waLBerla.tools',
            'waLBerla.tools.source_checker',
            'waLBerla.tools.report',
            'waLBerla.tools.sqlitedb',
            'waLBerla.tools.lbm_unitconversion',
            'waLBerla.tools.jobscripts']

collectFiles()

setup( name='waLBerla',
       version='1.0',
       author='Martin Bauer',
       author_email='martin.bauer@fau.de',
       url='http://www.walberla.net',
       packages=packages,
       package_data = {'' : ['walberla_cpp.so'] }
      )

import sys
if sys.argv[1] == 'build':
    print("\nCollected all files for waLBerla Python module.\n"
          " - to install run 'make pythonModuleInstall'\n"
          " - for development usage: \n"
          "      bash: export PYTHONPATH=%s:$PYTHONPATH \n"
          "      fish: set -x PYTHONPATH %s $PYTHONPATH \n" % ( walberla_binary_dir,walberla_binary_dir ) )