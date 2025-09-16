# Get and Build waLBerla  {#setup-instructions}

\tableofcontents

This guide describes the necessary steps to download and compile the waLBerla project sources.

## Prerequesites

To build waLBerla, you are going to need at least:

 - A recent C++20 compiler
 - CMake >= 3.24
 - An up-to-date MPI library for shared-memory parallel computing

WaLBerla also has numerous optional dependencies:

  - Nvidia CUDA for GPU computing on compatible GPUs
  - [Metis](https://github.com/KarypisLab/METIS) and [ParMetis](https://github.com/KarypisLab/ParMETIS)
    for static and dynamic load balancing
  - [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/) for handling of complex geometry meshes

> [!note]
> Many HPC clusters offer only outdated versions of CMake through their module system.
> For ways to install the latest CMake version in such a case, please take a look at
> [our FAQ](#faq-outdated-cmake).


## Download waLBerla's Source Code  {#get-sources}

### Clone the Repository

To work with the latest development revision of waLBerla,
clone the main [Git repository](https://i10git.cs.fau.de/walberla/walberla):

```bash
git clone https://i10git.cs.fau.de/walberla/walberla.git
```

If you mean to extend waLBerla's code base for your own purposes, or contribute to waLBerla's development,
you should create a fork off the main repository and clone the forked repo instead.

### Download a Release Distribution

To build a specific released version of waLBerla, you need to download its sources
as a Zip archive or tarball from [our GitLab](https://i10git.cs.fau.de/walberla/walberla/-/releases)
and extract that archive into your filesystem.

## Set Up the Build System

### Generate the Build System from a CMake Preset  {#cmake-preset}

WaLBerla's CMake build system exposes a large number of configuration options.
Also, depending on which of its modules you are using,
you may have to populate numerous additional build parameters with values specific to your local machine.
In order to document a build setup and to make builds reproducible,
we recommend to use [CMake Presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

Create a `CMakeUserPresets.json` file in the waLBerla project folder and populate it with the following
content:

<div class="tex2jax_ignore">

```json
{
  "version": 6,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 24,
    "patch": 0
  },
  "configurePresets": [
    {
      "name": "debug",
      "generator": "Ninja",
      "binaryDir": "${sourceDir}/build/${presetName}",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "Debug"
      },
      "environment": {
          "CC": "gcc",
          "CXX": "g++"
      }
    }
  ]
}
```

</div>

This defines a CMake *configure preset* called *debug*, using the `Ninja` generator
(swap this for `Unix Makefiles` if you don't have Ninja installed).
The CMake build type is set to `Debug`, and `gcc/g++` are used as compilers.

Generate the build system from this preset using the following command:

```bash
cmake --preset debug
```

This will set up the CMake build tree at the `build/debug` subfolder.

> [!note]
> IDEs like CLion or VS Code recognize CMake prefix files and permit you to
> select a preset for configuring their code analysis, auto-completion, and build tools.

> [!note]
> After modifying a preset you will have to regenerate the build system for changes to take effect.
> Either run `cmake --preset <preset-name>` again,
> or execute a regeneration task in your IDE
> (*Delete Cache and Reconfigure* in VS Code, *Reset Cache and Reload Project* in CLion).

### Validate the Build

To validate that the build system was correctly configured, navigate to the newly created subfolder `build/debug`.
There, run the following command to build the `UniformGridBenchmark` target to check if the framework successfully compiles:

```bash
cmake --build . --target UniformGridBenchmark
```

## Configure waLBerla's Build System

### Setting waLBerla Build Options

The build system of waLBerla is mostly configured through CMake cache variables.
If you are working on waLBerla directly, you should set these variables using the `cacheVariables` dictionary
of your CMake configure preset (see [above](#cmake-preset)), e.g.:

```json
"cacheVariables": {
  "CMAKE_BUILD_TYPE": "Debug",
  "WALBERLA_BUILD_WITH_OPENMP": true
}
```

Alternatively, cache variables can be set using the `-D` flag of the CMake command line.

If you are embedding waLBerla into another project, you can also use a preset or command line flags if
you wish to selectively enable certain features.
However, if your project needs to *always* enable or disable some feature, use the `option()` CMake function;
for instance, to always enable OpenMP in the waLBerla build, add the following to your CMakeLists:

```CMake
option( WALBERLA_BUILD_WITH_OPENMP "" ON )
```

## MPI

Support for MPI is enabled by default in waLBerla.
You may however disable it by setting the `WALBERLA_BUILD_WITH_MPI` build option to `false`.

### OpenMP

OpenMP is used in waLBerla for thread-based parallel computing on CPUs.
To enable it, set the `WALBERLA_BUILD_WITH_OPENMP` build option to `true`.

> [!warning]
> When enabling both OpenMP and MPI and running more than one MPI process per machine,
> take care not to oversubscribe your system.
> Use tools like [likwid-pin](https://github.com/RRZE-HPC/likwid/wiki/Likwid-Pin)
> or [likwid-mpirun](https://github.com/RRZE-HPC/likwid/wiki/Likwid-Mpirun) to control
> thread behavior in hybrid OpenMP/MPI applications.

## Build Tutorial, Benchmark and Showcase Apps

The waLBerla repository comes with a family of tutorials, benchmarks, and showcase applications, located in the `apps` directory.
These can be enabled in the build by setting the respective build options:

 - `WALBERLA_BUILD_TUTORIALS` for tutorials;
 - `WALBERLA_BUILD_BENCHMARKS` for benchmarks;
 - `WALBERLA_BUILD_SHOWCASES` for showcases.

Several tutorials, and most of the latest benchmark and showcase apps, require at least [code generation](#codegen-setup)
and the [embedded Python interpreter](#embedded-python) to be active in the build.
Also, [CUDA](#cuda-support) needs to be active for GPU applications.

## Automatic Code Generation  {#codegen-setup}

Most recent benchmark and showcase applications built with waLBerla are making heavy use of
automatic code generation using [pystencils](https://pycodegen.pages.i10git.cs.fau.de/pystencils/)
and [lbmpy](https://pycodegen.pages.i10git.cs.fau.de/lbmpy/).
Code generation takes place during the build process.

There are currently two separate code generation systems in place within waLBerla:
 - The legacy code generators, implemented in the `pystencils_walberla` and `lbmpy_walberla`
   Python modules, based on the pystencils and lbmpy 1.3 version line
 - The brand-new `SweepGen` code generator, based on pystencils 2.0 and meant to fully
   replace the former in the near future.

The two systems are implemented separately from one another and can coexists peacefully within a
build. However, you are encouraged to rely fully on `SweepGen` when building new applications,
as far as possible.

### SweepGen

SweepGen is the novel code generation system for numerical kernels in waLBerla.
To enable SweepGen in your build, set the `WALBERLA_ENABLE_SWEEPGEN` CMake option to `ON`;
e.g. in a preset:

```json
"cacheVariables": {
  "WALBERLA_ENABLE_SWEEPGEN": true
}
```

SweepGen requires an installation of Python >= 3.10.
Make sure that one is available, e.g. by setting the `Python_ROOT_DIR` variable to point
at the root of your Python installation (see also [FindPython](https://cmake.org/cmake/help/latest/module/FindPython.html)).

WaLBerla's build system will now automatically set up a virtual Python environment inside its build tree
and install all required packages.

For guidance on using SweepGen, refer to our [example apps](#example-apps) and the [SweepGen user manual](./sweepgen/index.html).

### Legacy Code Generation System

To enable the legacy code generators,
you need to manually install the required packages into your local Python environment
and point waLBerla's build system at the correct Python interpreter.

You are going to need at least Python 3.10.
If in doubt, run the query `python3 --version`
to find out the version of your Python installation.

#### Set Up the Python Environment

There are multiple providers for virtual Python environments out there.
We recommend using either the built-in [venv](https://docs.python.org/3/library/venv.html)
or the more versatile [virtualenv](https://virtualenv.pypa.io) module.
However, these only work if you already have an installation of Python 3.10 or newer on your machine;
which might not be the case on certain HPC clusters, or on outdated Linux distributions.
In that case, you could use [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)
to manage your Python environment.

##### Using venv or virtualenv

Create and activate a new [virtual environment](https://docs.python.org/3/library/venv.html)
inside your project directory (replace `-m venv` by `-m virtualenv` if using `virtualenv`):

```bash
python -m venv .venv
source .venv/bin/activate
```

Next, install pystencils and lbmpy (version 1.3), as well as jinja2 into that environment:

```bash
pip install pystencils~=1.3 lbmpy~=1.3 jinja2
```

##### Using Conda

Create a new Conda environment, using Python 3.10 and `pip`, located in the `conda-envs` subfolder of
your project directory. Then, activate the environment:

```
conda create -p ./conda-envs -n waLBerla-codegen python=3.10 pip
conda activate ./conda-envs/waLBerla-codegen
```

Then, use `pip` to install pystencils and lbmpy (version 1.3) and jinja2:

```bash
pip install pystencils~=1.3 lbmpy~=1.3 jinja2
```

> [!note]
> The pycodegen packages are not available on the conda installation channels,
> so you will have to install them using `pip`.
> For an overview of pitfalls when using `pip` inside conda environment, refer to
> [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#using-pip-in-an-environment).

#### Configure CMake for Code Generation

There are two cache variables you need to set in order to activate build-time code generation
in CMake:

 - Set `WALBERLA_BUILD_WITH_CODEGEN` to `true`
 - Set `Python_ROOT_DIR` to the `bin` folder of your Python environment.

When using `venv` or `virtualenv` according to the above instructions, `Python_ROOT_DIR`
will be `<project-dir>/.venv/bin`.
When using `conda` like above, it will be `<project-dir>/conda-envs/waLBerla-codegen/bin`.

In the presets file, you can set these options like this:

<div class="tex2jax_ignore">

```json
"cacheVariables": {
  "WALBERLA_BUILD_WITH_CODEGEN": true,
  "Python_ROOT_DIR": "<path-to-python-root>"
}
```

</div>

Here is a minimal working CMake preset for code generation, using `venv`:

<div class="tex2jax_ignore">

```json
{
  "name": "minimal-codegen",
  "generator": "Ninja",
  "binaryDir": "${sourceDir}/build/${presetName}",
  "cacheVariables": {
    "WALBERLA_BUILD_WITH_CODEGEN": true,
    "Python_ROOT_DIR": "${sourceDir}/.venv/bin"
  }
}
```

</div>

## Embedded Python Interpreter  {#embedded-python}

WaLBerla offers the option to run Python scripts in an embedded Python interpreter to configure and
parametrize simulations.
This requires that both CPython >= 3.10 and the Python development headers are installed.
When this is the case, you may activate the Python coupling by setting the `WALBERLA_BUILD_WITH_PYTHON`
build option to `true`.

To use the embedded Python interpreter, waLBerla needs the Python package pybind11, which is usually fetched automatically.

But there are systems, on which pybind11 can not be fetched automatically, because no network connection is possible.
In this case, download the version 2.13.6 of the pybind11 source from https://github.com/pybind/pybind11/tree/v2.13,
copy it to the system without internet connection and set the CMake option ’CMAKE_PYBIND11_SOURCE_DIR’ to your local pybind11 source directory.


## CUDA Support {#cuda-support}

### Prerequisites

For a CUDA support in CMake, installing the NVIDIA GPU drivers and the CUDA Toolkit is required.
You can verify your installation by running

```bash
# test if CUDA compiler is accessible
which nvcc
nvcc --version

# monitor NVIDIA GPU and check if drivers are set up
nvidia-smi
```

<div class="tex2jax_ignore">

> [!note]
> In case `nvcc` is not found, adding it to `PATH` may resolve the issue:
>
> ```
> # may be located somewhere else depending on installation
> export CUDA_HOME=/usr/local/cuda
>
> export PATH=$PATH:$CUDA_HOME/bin
> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib64
> ```

</div>

### Configure CMake for CUDA

The CUDA support in CMake is enabled via `WALBERLA_BUILD_WITH_CUDA`.
When this flag is set to `ON`, waLBerla automates several steps for using CUDA in the CMake project.
It attempts to find and configure CUDA, enables compile options, sets up preprocessor macros such as `WALBERLA_BUILD_WITH_GPU_SUPPORT`, etc.

It also configures the build pipeline such that the CUDA sources are compiled separately with the compiler specified in `CMAKE_CUDA_COMPILER`.

The `CMAKE_CUDA_ARCHITECTURES` specifies the compute capability of the target GPU architecture,
which is then used by `nvcc` to generate hardware-specific code.
If you are not sure which compute capability your GPU has,
refer to [the NVidia Developer Documentation](https://developer.nvidia.com/cuda-gpus).

> [!note]
> Compute capability X.Y translates to `CMAKE_CUDA_ARCHITECTURES=XY`, e.g. 8.0 means `CMAKE_CUDA_ARCHITECTURES=80`.

You can also use `CMAKE_CUDA_ARCHITECTURES=native` to have the compute capability determined
automatically, based on your local GPU.

To use CUDA with waLBerla, set at least `WALBERLA_BUILD_WITH_CUDA` and `CMAKE_CUDA_ARCHITECTURES`, e.g. in the CMake presets file:

<div class="tex2jax_ignore">

```json
"cacheVariables": {
  "CMAKE_CUDA_COMPILER": "$env{CUDA_HOME}/bin/nvcc",
  "WALBERLA_BUILD_WITH_CUDA": true,
  "CMAKE_CUDA_ARCHITECTURES": "61",
}
```

</div>

Or on the command line:

```
cmake \
    -DCMAKE_CUDA_COMPILER=$CUDA_HOME/bin/nvcc \
    -DWALBERLA_BUILD_WITH_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=61 \
    ...
```

## OpenMesh for Geometry Processing

WaLBerla's `mesh_common` and `mesh` modules couple with [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download),
an open-source tool developed at RWTH Aachen,
to import and handle surface meshes of geometric objects.

The OpenMesh tool is fetched automatically by CMake, when setting `WALBERLA_BUILD_WITH_OPENMESH=ON`.

There are systems, on which OpenMesh can not be fetched automatically, because no network connection is possible.
In this case, download the newest OpenMesh source from https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git, 
copy it to the system without internet connection and set the CMake option option ’CMAKE_OPEN_SOURCE_DIR’ to your local OpenMesh source directory.

## Metis and ParMetis for Load Balancing

WaLBerla optionally links against the popular graph partitioning libraries [Metis](https://github.com/KarypisLab/METIS)
and [ParMetis](https://github.com/KarypisLab/ParMETIS).
To activate Metis and ParMetis support, set the respective build options to `true`:

 - `WALBERLA_BUILD_WITH_METIS` for Metis
 - `WALBERLA_BUILD_WITH_PARMETIS` for ParMetis

The build system will then attempt to find a local installation of either library.

> [!note]
> To use (Par)Metis with waLBerla, Metis must have been built with 64-bit integers and floats enabled
> (see [Metis Build Configuration Options](https://github.com/KarypisLab/METIS?tab=readme-ov-file#common-configuration-options-are)).


