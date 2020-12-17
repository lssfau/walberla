import sys
import platform
from os.path import exists, join
import shutil

from setuptools import setup

# The following variables are configure by CMake
walberla_source_dir = "${walberla_SOURCE_DIR}"
walberla_binary_dir = "${CMAKE_CURRENT_BINARY_DIR}"
suffix = "${PYTHON_MODULE_EXTENSION}"
prefix = "${PYTHON_MODULE_PREFIX}"
walberla_module_file_name = prefix + "walberla_cpp" + suffix

version = "${WALBERLA_VERSION}"

if platform.system() == 'Windows':
    configuration = 'Release'
else:
    configuration = ''


def collectFiles():
    src_shared_lib = join(walberla_binary_dir, configuration, walberla_module_file_name)
    dst_shared_lib = join(walberla_binary_dir, 'waLBerla', walberla_module_file_name)
    # copy everything inplace

    if not exists(src_shared_lib):
        print("Python Module was not built yet - run 'make walberla_cpp'")
        exit(1)

    shutil.rmtree(join(walberla_binary_dir, 'waLBerla'), ignore_errors=True)

    shutil.copytree(join(walberla_source_dir, 'python', 'waLBerla'),
                    join(walberla_binary_dir, 'waLBerla'))

    shutil.copy(src_shared_lib,
                dst_shared_lib)


packages = ['waLBerla',
            'waLBerla.evaluation',
            'waLBerla.tools',
            'waLBerla.tools.source_checker',
            'waLBerla.tools.report',
            'waLBerla.tools.sqlitedb',
            'waLBerla.tools.lbm_unitconversion',
            'waLBerla.tools.jobscripts',
            'waLBerla.tools.config']

collectFiles()


setup(name='waLBerla',
      version=str(version),
      author='Markus Holzer',
      author_email='markus.holzer@fau.de',
      url='http://www.walberla.net',
      packages=packages,
      package_data={'': [walberla_module_file_name]}
      )

if sys.argv[1] == 'build':
    print("\nCollected all files for waLBerla Python module.\n"
          " - to install run 'make pythonModuleInstall'\n"
          " - for development usage: \n"
          "      bash: export PYTHONPATH=%s:$PYTHONPATH \n"
          "      fish: set -x PYTHONPATH %s $PYTHONPATH \n" % (walberla_binary_dir, walberla_binary_dir))
