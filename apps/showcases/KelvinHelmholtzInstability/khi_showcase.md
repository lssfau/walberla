# Kelvin-Helmholtz Instability App – Showcase
---

## Overview

This app simulates a quasi-2D single-mode Kelvin-Helmholtz instability for a range of LBM schemes, Reynolds numbers and grid resolutions.  The primary quantitative quantities of interest are the momentum thickness and average kinetic energy.  Qualitatively, isosurfaces of a passive scalar are used to visualise  vorticity-carrying eddies, which provides further insight into the flow symmetry. This showcase is therefore an example of multi-distribution one-way coupling between a hydrodynamic solver and an advection-diffusion solver.

---

## 1. Configuring the App
By default, the app is configured to run without requiring any user input. The CMakeLists.txt file iterates over the LBM configurations of interest, invoking the KHI.py code-generation script for each configuration to generate the corresponding binaries.

To add new configurations, follow these steps:
1. Update CMakeLists.txt: Add a new config ID (e.g., test_config) to the relevant section of the file.
2. Modify KHI.py: Extend the options_dict in KHI.py with the details for the new LBM configuration. For example:

To add additional cases, one must update the `CMakeLists.txt` with an appropriate config ID tag (i.e., `test_config`), then create associated LBM configuration details in the `KHI.py`  `options_dict`. For example:

```
options_dict = {
    'test_config': {
        'subgrid_scale_model': ....,
        'fourth_order_correction': ....,
    },
    ...
}
```
Once these changes are made, simply rebuild the app and run it!

## 2. Building the App
### General Build
The cmake flags required for build this showcase (on CPU) are:
- Set `WALBERLA_BUILD_SHOWCASES` to `ON`
- Set `WALBERLA_BUILD_WITH_CODEGEN` to `ON`
- Set `WALBERLA_BUILD_WITH_PYTHON` to `ON`

### GPU Build Instructions
For a NVIDIA GPU, the additional cmake flags must be specified:
- Set `WALBERLA_BUILD_WITH_CUDA` to `ON`
- Set `CMAKE_CUDA_ARCHITECTURES` (GPU Specific)

For a AMD GPU, the additional cmake flags must be specified:
- Set `WALBERLA_BUILD_WITH_HIP` to `ON`

### Example CMake Configuration
For a NVIDIA GeForce RTX 4090 
```
cmake \
    -DWALBERLA_BUILD_SHOWCASES=ON \
    -DWALBERLA_BUILD_WITH_CODEGEN=ON \
    -DWALBERLA_BUILD_WITH_PYTHON=ON \
    -DCMAKE_CUDA_COMPILER=$CUDA_HOME/bin/nvcc \
    -DWALBERLA_BUILD_WITH_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=89 \
    /path/to/walberla/
```

### Compilation
```bash
cd apps/showcases/KelvinHelmholtzInstability/ && make -j
```

## 3. Parameterising and Running the App
### Parameterisation
Acoustic scaling is automatically applied within the KHI.cpp app, converting all SI quantities defined in `KHI_parameters.py` into their lattice equivalent. All physical quantities are denoted by the `_SI` suffix.

All modifications to the case setup or physical parameters should be implemented in `KHI_parameters.py`.

## Running the App
Single-Core / Single-GPU Execution:
```bash
./KHI_<config> KHI_parameters.py
```
Multi-Core / Multi-GPU (MPI) Execution:
```bash
mpirun -np <nProc> ./KHI_<config> KHI_parameters.py
```
Replace `<config>` with the CMakeList config identifier and `<nProc>` with the desired number of processes.

## 4. Validation
Due to the highly non-linear nature of this problem, no analytical solution exists. Consequently, validation can be performed by:

- Computing a DNS (Direct Numerical Simulation) reference solution and comparing the app's results against this high-fidelity baseline.
- Benchmarking against other solvers or published results available in the literature.

Common metrics for evaluating the performance are the non-dimensional momentum thickness and plane-averaged kinetic energy:

**Quantitative (provided as outputs in the code):**

- **Momentum thickness:**

  $\Theta^* = \frac{\Theta}{L}$

- **Plane-averaged kinetic energy:**

  $k^* = \frac{\overline{u_i u_i}}{2 \overline{U}^2}$

These quantities are often plotted against the vortex timescale:

$\tau^* = \frac{t \Delta U}{L}$

Further qualitative insight can be achieved through analysing the isosurface contour of the passive scalar.

## 5. General Notes
- Gaussian quadrature is valid for smooth functions.  which is the case when the initial contact surface lies exactly halfway between two adjacent cell centers. This condition is met when an even number of cells is used in each direction.  Using an odd number of cells on a coarse grid may lead to spurious results.  On refined grids the initial perturbation growth can be suppressed.
- The representation of the hyperbolic-tangent mean-velocity profile on highly-refined grids is smooth. This may lead to an increased mometum thickness at the early timescales, or in some cases, prevent the Kelvin-Helmholtz instability from triggering altogether. To mitigate this, set a smaller value for `momentumThickness_SI` in `KHI_parameters.py`.

## 6. References and Further Reading
TODO: Once published
- This app is based on the configuration used in doi:XXX
- Other useful references include:
    - doi:10.1002/fld.3853
    - doi:10.1093/mnras/stv2564








