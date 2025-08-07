# waLBerla

waLBerla (widely applicable Lattice Boltzmann from Erlangen) is a massively parallel framework for multiphysics simulation applications.
Beyond computational fluid dynamics with the lattice Boltzmann method, the framework now features multiphase and free-surface flows,
rigid body and particle dynamics as well as fluid-structure coupling with moving geometries.
It scales from laptops to current and future supercomputers while maintaining
near-perfect efficiency.

## Getting Started with waLBerla

Refer to our [Setup Guide](https://walberla.pages.i10git.cs.fau.de/walberla/setup-instructions.html) for instructions
on setting up and building waLBerla.

## Documentation and Tutorials

You can find our framework documentation, guides, tutorials, and examples on the following pages:

 - Latest Release: [C++ Framework](http://walberla.net/doxygen/index.html), [Python Interface](http://walberla.net/sphinx/index.html)
 - Current Development Revision: [C++ Framework](walberla.pages.i10git.cs.fau.de/walberla/)

## Get involved

### Contributing

Please refer to the [contribution guide](CONTRIBUTING.md) for guidance on contributing to waLBerla.

### Support

To get in touch with the waLBerla developers, use our [Issue Tracker](https://i10git.cs.fau.de/walberla/walberla/issues)
or the waLBerla mailing list ([cs10-walberla@fau.de](mailto:cs10-walberla@fau.de)).

## Authors

Many thanks go to waLBerla's [contributors](AUTHORS.txt)

### Please cite us

If you use waLBerla in a publication, please cite the following articles:

Overview:
- M. Bauer et al., *waLBerla: A block-structured high-performance framework for
  multiphysics simulations*. Computers & Mathematics with Applications, 2020.
  https://doi.org/10.1016/j.camwa.2020.01.007.

Grid Refinement:
- F. Schornbaum and U. Rüde, *Massively parallel algorithms for the lattice Boltzmann
  method on nonuniform grids*. SIAM Journal on Scientific Computing, 2016.
  https://doi.org/10.1137/15M1035240

LBM - Particle Coupling:
- C. Rettinger and U. Rüde, *A comparative study of fluid-particle coupling methods for
  fully resolved lattice Boltzmann simulations*. Computers & Fluids, 2017.
  https://doi.org/10.1016/j.compfluid.2017.05.033

Free-surface LBM:
- C. Schwarzmeier et al., *Comparison of free-surface and conservative Allen-Cahn phase-field
  lattice Boltzmann method*. Journal of Computational Physics, 2023.
  https://doi.org/10.1016/j.jcp.2022.111753

Allen-Cahn phase-field LBM
- M. Holzer et al., *Highly efficient lattice Boltzmann multiphase simulations of immiscible
  fluids at high-density ratios on CPUs and GPUs through code generation*. The International Journal of High Performance Computing Applications, 2021.
  https://doi.org/10.1177/10943420211016525

MESA-PD:
- S. Eibl and U. Rüde, *A Modular and Extensible Software Architecture for Particle Dynamics*.
  Proceedings of the 8th International Conference on Discrete Element Methods.
  https://mercurylab.co.uk/dem8/full-papers/#page-content

Carbon Nanotubes:
- G. Drozdov et al., *Densification of single-walled carbon nanotube films:
  Mesoscopic distinct element method simulations and experimental validation*.
  Journal of Applied Physics, 2020. https://doi.org/10.1063/5.0025505

## License

waLBerla is licensed under [GPLv3](COPYING.txt).
