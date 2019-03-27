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

The minimum requirements are a C++14-compliant compiler (e.g. GCC or Clang),
the [Boost](http://www.boost.org) library and the [CMake](http://www.cmake.org)
build system. Furthermore, you need an MPI library (like
[Open MPI](http://www.open-mpi.org)) if you want to make use of parallel
processing capabilities. All of these dependencies are typically available in
your operating system's package manager.

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

If you use waLBerla in a publication, please cite the following article:

- C. Godenschwager, F. Schornbaum, M. Bauer, H. Köstler, and U. Rüde. A
framework for hybrid parallel flow simulations with a trillion cells in complex
geometries. In: Proceedings of the International Conference on High Performance
Computing, Networking, Storage and Analysis, page 35. ACM, 2013.

## License

waLBerla is licensed under [GPLv3](COPYING.txt).
