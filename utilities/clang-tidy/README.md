# Using clang-tidy on waLBerla {#clang-tidy}

See also the [Clang-Tidy Documentation](https://clang.llvm.org/extra/clang-tidy/).

This document briefly explains how to run a clang-tidy analysis and (auto-)fixing pass on the waLBerla code base.

## tl;dr

 1. Set up build system: `cmake --preset clang-tidy`
 2. Navigate to `build/clang-tidy`
 3. Run this command: `python utilities/clang-tidy/analyze.py -p utilities/clang-tidy/analyze.yml -c compile_commands.json -r ../..`
 4. Inspect the results
 5. Fix the errors; if applicable, use the autofixing feature.


## Set up CMake

To run clang-tidy, CMake needs to be configured to emit a compile command data base which clang-tidy
reads to figure out which files, with which compilation options, it should analyze.
Also, the build system should be set up in Debug mode with additional debugging features enabled
such that all debug code paths can be covered.

Here is a sample `CMakeUserPresets.json` with a configure preset called `clang-tidy`, which defines
a build system that can be used for clang-tidy analysis:

```JSON
{
  "version": 6,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23,
    "patch": 0
  },
  "configurePresets": [
    {
      "name": "clang-tidy",
      "generator": "Unix Makefiles",
      "binaryDir": "${sourceDir}/build/clang-tidy",
      "cacheVariables": {
        "CMAKE_EXPORT_COMPILE_COMMANDS": true,
        "WALBERLA_BUFFER_DEBUG": true,
        "WALBERLA_BUILD_TESTS": true,
        "WALBERLA_BUILD_BENCHMARKS": true,
        "WALBERLA_BUILD_TUTORIALS": true,
        "WALBERLA_BUILD_TOOLS": true,
        "WALBERLA_BUILD_WITH_MPI": true,
        "WALBERLA_BUILD_WITH_OPENMP": true,
        "CMAKE_BUILD_TYPE": "Debug",
        "WALBERLA_BUILD_WITH_METIS": true,
        "WALBERLA_BUILD_WITH_PARMETIS": true,
        "WALBERLA_BUILD_WITH_OPENMESH": true,
        "WALBERLA_DOUBLE_ACCURACY": true,
        "WALBERLA_LOGLEVEL": "DETAIL"
      }
    },
  ]
}
```

The above configuration:
 - requires that *OpenMesh* is installed - if you do not plan on analyzing the `mesh_common` and `mesh` modules,
   you can set `WALBERLA_BUILD_WITH_OPENMESH` to `false`.
 - requires that *Metis* and *ParMetis* are installed - if you do not plan to analyze the Metis-dependent code
   in `core` and `blockforest`, you can set `WALBERLA_BUILD_WITH_METIS` and `WALBERLA_BUILD_WITH_PARMETIS` to `false`.
 - prepares the build system also for tests, benchmarks and tutorials - if you do not wish to analyze those, you can exlude them.

This preset is contained in the `CMakePresets.json` file. To generate the build system defined by it, run `cmake --preset clang-tidy`.

## Run the Analysis

We provide the `analyze.py` script in this directory to automate the clang-tidy analysis of the waLBerla code base.
For an overview of its options, run `python analyze.py --help`.

To run the full analysis pass, with your shell located in the waLBerla project root directory, you can use the following command:

```
python utilities/clang-tidy/analyze.py -p utilities/clang-tidy/analyze.yml -c build/clang-tidy/compile_commands.json
```

In there, the `-p` option specifies the parameter file, and the `-c` option says where to find the compile database.
You can specify a custom parameter file to restrict the set of modules and apps to be analyzed.
