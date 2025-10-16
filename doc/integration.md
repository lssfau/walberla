# Using waLBerla as a library in your CMake Project  {#using-as-library}

The recommended way to embed waLBerla into your own CMake project is to integrate its build system
with your own CMake tree, such that waLBerla is built from source alongside your project.
There are two ways to achieve this in a stable and reproducible manner: [Git submodules](#git-submodule)
and [FetchContent](#fetchcontent)

## Include waLBerla as a Git Submodule  {#git-submodule}

To include waLBerla as a submodule, run the following commands inside your project folder:

```bash
git submodule add https://i10git.cs.fau.de/walberla/walberla.git
git submodule update --init
```

After waLBerla has been downloaded, `git commit ` the submodule.
You may then add waLBerla to your CMake project by adding the following to your root `CMakeLists.txt`:

```CMake
add_subdirectory( walberla )
```

## Pull waLBerla using FetchContent  {#fetchcontent}

CMake allows us to download and include external projects into a build tree at configuration time using `FetchContent`.
This can be more convenient and easier to reproduce than Git submodules, or other ways of semi-manual installation.

To include waLBerla into your CMake project using `FetchContent`, add the following lines to your `CMakeLists.txt`:

```CMake
include(FetchContent)

FetchContent_Declare(
    walberla
    GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla.git
    GIT_TAG <commit-hash / branch name / git tag>
)

FetchContent_MakeAvailable(walberla)
```

Using the `GIT_TAG` option you can specify the exact Git revision to include:
We recommend fixing this to a specific commit hash or tag to avoid surprising changes.
In this case, you will have to remember to regularily update this dependency.

## Manage waLBerla's build options

How to manage waLBerla's CMake options in a consumer project depends on whether or not a
certain option is allowed to change in the context of your project, or if its value should be
fixed.
In the former case, set the option through `-D` CMake command line arguments or in a preset as usual.
In the latter case, override the option's value by using the `option` command *before* including waLBerla;
e.g. to always enable OpenMP:

```CMake
option( WALBERLA_BUILD_WITH_OPENMP "" ON )

# ...

add_subdirectory( walberla )
```
