# waLBerla

waLBerla (widely applicable Lattice Boltzmann from Erlangen) is a massively
parallel framework for multi physics applications. Besides its original
objective, Lattice Boltzmann solvers for hydrodynamics, it now contains
modules for other applications like Multigrid and rigid body dynamics
as well. Great emphasis is placed on the interoperability between the modules
in particular the fluid-particle coupling.
It scales from laptops to current and future supercomputers while maintaining
near-perfect efficiency.

See https://www.walberla.net/ for more information and a showcase of applications.

## Documentation and Tutorials

Documentation for the C++ framework is available in
[Doxygen](http://walberla.net/doxygen/index.html), while the Python interface
is documented in [Sphinx](http://walberla.net/sphinx/index.html).

## Getting started

The minimum requirements are a C++17-compliant compiler (e.g. GCC or Clang)
and the [CMake](http://www.cmake.org)
build system. Furthermore, you need an MPI library (like
[Open MPI](http://www.open-mpi.org)) if you want to make use of parallel
processing capabilities. All of these dependencies are typically available in
your operating system's package manager.

### CMake

The typical steps, assuming your are in the waLBerla source directory, are:

- `mkdir build; cd build` create a build directory and change into it
- `cmake ..` call CMake with the waLBerla source directory as an argument
- `make` build waLBerla

To specify a CMake option you need to use `-D(Option)=(Value)`. For example to set the C++ compiler one can use:
`cmake -DCMAKE_CXX_COMILER=clang++`

To list and modify the CMake options the `ccmake` tool can be used. Just call `ccmake .` in your **build** directory to see and change the
CMake options and variables.

Some important CMake variables:

- `WALBERLA_BUILD_WITH_CODEGEN` Enable pystencils code generation"
- `Python_ROOT_DIR` Specify the directory of the `python` executable. e.g. `~/miniconda/bin/`
- `MPI_HOME` Specify the base directory of the MPI installation.
- `WALBERLA_BUILD_WITH_PYTHON` Support for embedding Python
- `WALBERLA_BUILD_WITH_CUDA` Enable CUDA support

For a full list of CMake Option see the [CMakeLists.txt](CMakeLists.txt) file or use `ccmake` as described above.

### Codegen and Python

To use the `lbmpy`/`pystencils` code generation please install the packages with e.g. `pip3 install lbmpy` and specify the correct python
environment when calling CMake.

In previous versions of CMake one could use `PYTHON_EXECUTABLE` or `PYTHON_ROOT_DIR` (all upper case) to specify the python executable or
the directory. This does **NOT** work anymore. Please use `Python_ROOT_DIR`.

## Get involved

### Contributing

Please submit all code contributions on our
[GitLab](https://i10git.cs.fau.de/walberla/walberla). To get an account, please
sign and submit the [contributor license agreement](CONTRIBUTING.txt).

### Support

While we currently do not have a mailing list, any questions can be asked via
the [Issue Tracker](https://i10git.cs.fau.de/walberla/walberla/issues).

## Authors

Many thanks go to waLBerla's [contributors](AUTHORS.txt)

### Please cite us

If you use waLBerla in a publication, please cite the following articles:

Overview:
- M. Bauer et al, *waLBerla: A block-structured high-performance framework for
  multiphysics simulations*. Computers & Mathematics with Applications, 2020.
  https://doi.org/10.1016/j.camwa.2020.01.007.

Grid Refinement:
- F. Schornbaum and U. Rüde, *Massively parallel algorithms for the lattice boltzmann
  method on nonuniform grids*. SIAM Journal on Scientific Computing, 2016.
  https://doi.org/10.1137/15M1035240

LBM - Particle Coupling:
- C. Rettinger and U. Rüde, *A comparative study of fluid-particle coupling methods for
  fully resolved lattice Boltzmann simulations*. Computers & Fluids, 2017.
  https://doi.org/10.1016/j.compfluid.2017.05.033

MESA-PD:
- S. Eibl and U. Rüde, *A Modular and Extensible Software Architecture for Particle Dynamics*.
  Proceedings Of The 8Th International Conference On Discrete Element Methods.
  https://mercurylab.co.uk/dem8/full-papers/#page-content

Carbon Nanotubes:
- G. Drozdov et al, *Densification of single-walled carbon nanotube films:
  Mesoscopic distinct element method simulations and experimental validation*.
  Journal of Applied Physics, 2020. https://doi.org/10.1063/5.0025505

## License

waLBerla is licensed under [GPLv3](COPYING.txt).
