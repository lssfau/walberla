# SweepGen User Manual

Welcome to the waLBerla SweepGen User Manual.

SweepGen is the novel metaprogramming system for numerical kernels in waLBerla.
It integrates the code generation capabilities of [pystencils][pystencils]
and the mathematical modelling framework of [lbmpy][lbmpy]
with waLBerla's data structures and algorithmic idioms.
SweepGen is a Python library that is operated through user-defined *generator scripts*
written in Python.
Its primary purpose is to generate the eponymous *Sweep* classes
which encompass parallel numerical update rules,
but it also produces a number of supporting components,
such as *pack-info* classes for ghost-layer exchange.

## Quickstart

To activate SweepGen for your waLBerla build,
set the `WALBERLA_ENABLE_SWEEPGEN` CMake variable to true, e.g.:

```{code-block} bash
cmake -S . -B build -DWALBERLA_ENABLE_SWEEPGEN=ON
```

The build system will now automatically set up a virtual Python environment
and install all packages required for SweepGen.

SweepGen runs through generator scripts during the build process of your waLBerla app.

## Contents

:::{toctree}
:caption: How-To Guides
:maxdepth: 1

guides/IDEandDebugging
:::

:::{toctree}
:caption: Reference
:maxdepth: 1

reference/PythonEnvironment
reference/Configuration
reference/SweepGeneration
reference/GeometryFunctions
reference/Prefabs
:::


[pystencils]: https://pycodegen.pages.i10git.cs.fau.de/docs/pystencils/2.0dev/
[lbmpy]: https://pycodegen.pages.i10git.cs.fau.de/lbmpy/
