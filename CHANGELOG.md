# Changelog

## [6.1] - 2022-07-25
### Added
- Free-surface LBM extension:
  - Add implementation
  - Add several showcases
  - Add several tests
- LBM - MESA_PD coupling:
  - Add partially saturated cells method (PSM)
  - Add fluidized bed showcase
  - Add virtual mass stabilization technique for light particles
  - Add support for more shapes, e.g., convex polyhedron
- MESA_PD:
   - Add extensive application for dense particle packing generation
- AMD - HIP support
  - Support of the ROCm Toolchain and thus AMD HIP as second GPU language
  - All CUDA related files, namespaces, folders etc are renamed to gpu.
  - Include "GPUWrapper.h" to use general GPU functions cudaMalloc -> gpuMalloc
  - WALBERLA_BUILD_WITH_HIP and WALBERLA_BUILD_WITH_GPU_SUPPORT as new CMake variables introduced

### Changed
- Update and extend phase-field LBM showcases
- Allow access to PDF centering information (for being used in generated LBM kernels)
- Adapt code generation backend to be compatible with pystencils 1.0 and lbmpy 1.0
- Required minimum dependencies:
  - C++17-compliant compiler
  - CMake 3.14
  - pybind 2.6.2
  - lbmpy 1.0
  - pystencils 1.0

### Deprecated
- GUI

## [5.1] - 2020-04-09
### Added
- Add new tutorials and showcases
- Extend MESA-PD functionalities, including several molecular dynamics models
- Fluid-particle coupling with MESA-PD: functionalities, tests, benchmark scenarios

### Changed
- Update to C++17
- Update CUDA compiler support
- Extend Clang-Tidy coverage
- Add closer integration of code generation using pystencils and lbmpy
- Python Coupling now build upon pybind11. Boost.Python is no longer supported
  - lbm module dropped from python coupling due to deprecation for a long time
  - geometry, postprocessing and timeloop dropped from python coupling due to its low usage
  - PEP8-ification of Python API. This means all keyword arguments are now in snake_case and not in CamelCase as before.

### Fixed
- Guo force model for non-SRT, may change simulation results

## [4.1] - 2019-04-19
### Added
- Galerkin coarsening for Multigrid
- LBM-PE-Coupling:
  - new coupling approach for unresolved particle interactions (discrete particle method)
  - adaptive grid refinement for coupled simulations
  - load balancing functionalities for coupled simulations
  - module description
- integrated *pystencils* and *lbmpy* code generation for kernels and pack infos
- new GPU communication, including support for GPUDirect
- load balancing functionality for the pe
- implemented IProbe communication as an alternative to two message communication for unknown size communication
- new creation helpers for BlockForest
- Minor:
   - dynamic load balancing now possible with levels ignored
   - `ExtendedBoundaryHandlingFactory` now uses `ParserUBB` instead of `UBB` so that velocity profiles can be specified as an equation in a parameter file
   - Enabled the body selection functions for PSM coupling method
   - grid_generators now allow range based for loops

### Changed
- A compiler with full C++14 support is now required
- All Boost usage has been replaced with the corresponding standard library functionality, except for Boost.Python (used for the `python_coupling` module), Boost.PropertyTree (used in `config::configToBoostPropertyTree`) and Boost.Graph (used by `math::EquationSystem`). This usually means you need to replace `boost::` with `std::` in your code and change some `#include`s.
- API changes in blockforest::PhantomBlockForest, blockforest::loadbalancing, pe::amr::weight_assignment
- API change for vtk::VTKOutput::getFilenames
- made SendBuffer memory access more high level
- PE coupling:
   - changed body mapping functions: removed most of them, added mapping-decider function, accomodated changes to test cases, added new mapping test case
   - changed pe coupling mapping functions interfaces, added new mapping functions, adapted test cases and benchmarks
   - change in lubrication correction functionality to not require the lattice model but use the dynamic viscosity directly
- PE:
   - rebased Union on BodyStorage
   - `pe::Union`, `boundary::BoundaryHandling` and `boundary::BoundaryHandlingCollection` are now variadic templates instead of taking a tuple of bodies/boundaries/handlers. This means that you need to replace `std::tuple<A,B>` with `A,B` in these cases.
   - using smart pointers for all memory management
   - made setMass protected
   - extended DEM collision model with dt
   - pe::createBlockForest changed to support a variable number of processes
   - changed BodyStatistics to use shared_ptr instead of reference

### Removed
- Remove dependency resolution from Singleton
- from PE
   - Node
   - PtrVector
   - contacts from RigidBody
   - attached bodies
   - attachables

### Fixed
- Fix implict/explict naming in pe::cr

### Deprecated
- all dynamic level-wise balance functions (use the more general ones, without "level-wise")
