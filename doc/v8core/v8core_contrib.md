# V8 Core Library Contributor's Guide {#v8core-contrib}

Welcome to the contributors guide to the waLBerla v8 core library.
Here you will find instructions and guidance on how to contribute code to the
v8 core modules.
Please read this guide carefully before beginning development.

## Tenets

The *V8 core library* is part of the waLBerla project's ongoing efforts at modernization.
It is designed as a platform for the *reimplementation* of waLBerla's core framework features,
and will replace the legacy framework modules over time.

In developing the V8 core, we follow these guiding principles:
 - *Portability First:* The V8 core library is designed with performance-portability as a fundamental goal.
   We accomodate both CPU and GPU targets as first-class citizens,
   and design API abstractions that permit the transparent usage of accelerators with minimal code overhead,
   while building the underlying logic to be as efficient as possible.
 - *Bleeding Edge Coding:* We make use of state-of-the-art C++ coding practice, using features of the C++20
   standard (such concepts, ranges, `constexpr` programming, etc) to maximum effectivenes.
   In doing so, we accept incompatibility with older compiler versions.
 - *Modular and Layered Architecture:* We design a strictly modular architecture, following the single-responsibility principle.
   We prefer small, self-contained components with clear interfaces over larger, monolithic structures.
   We identify suitable architecture layers, providing interfaces at various layers of abstraction to accomodate
   differing levels of user requirements and expertise.
 - *Rigorous Quality Assurance:* We verify code correctness and stability with extensive unit
   and integration tests, and apply static code analysis to ensure a high level of code quality.

## Structure

### Source Files

The library's public header files are located in the directory `include/walberla/v8`,
placed in subfolders according to module affiliation.
We employ the file extension `.hpp` and `#pragma once` include guards for headers.

Translation units are located at the `src/v8` directory.

Each source file must have a license header, including authorship information.
Source files must be formatted using `clang-format` according to the style laid down in
the waLBerla project's `.clang-format` file.

### Modules and Namespaces

The outermost namespace is `walberla::v8`.
Each module `ModuleX` has its own namespace `walberla::v8::module_x`.
Module namespace names are written in lower snake-case.
A module's components and submodules should be further distributed into child namespaces, as is appropriate.
Modules are allowed to export their most important APIs to the parent namespace `walberla::v8`
via `using`-declarations.

The header files and translation units of module `ModuleX` are located in the `v8/module_x` subfolders
of the `include/walberla` and `src` directories, respectively.
Each module has a primary header file `walberla/v8/ModuleX.hpp` that `#include`s and re-exports
the module's public APIs.
The primary header files are in turn included the bulk-header `walberla/V8.hpp`.

### Documentation

Public APIs must be sufficiently documented.
We use [Doxygen documentation comments](https://www.doxygen.nl/manual/docblocks.html) in Javadoc-style
for API documentation.
We use the `@` syntax for special commands, and prefer Markdown over Doxygen-commands and raw HTML
(see [Markdown in the Doxygen documentation](https://www.doxygen.nl/manual/markdown.html)).

Classes and free functions must be *grouped* according to their owning modules and submodules,
in order to be easily locatable in the documentation pages.
See [Grouping in the Doxygen documentation](https://www.doxygen.nl/manual/grouping.html#topics).

 - In `walberla/V8.hpp`, the master group is defined:
   ```
   @defgroup v8core
   ```
 - The primary header file `ModuleX.hpp` of each top-level module `ModuleX`
   must define the subgroup `v8core-moduleX` like this:
   ```
   @defgroup v8core-moduleX Module Title
   @ingroup v8core
   @brief Module Short Description

   ... Module Long Description ...

   ```
 - Module documentation should be split further into sub-groups according to their components.

When done correctly, API documentation will appear cleanly structured in the [API Reference > Modules](#v8core) section.

### Testing {#v8-contrib-testing}

Code in the V8 core modules must be rigorously tested.
In waLBerla, unit and integration testing is orchestrated using [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html).
For writing tests, we have developed a [toolkit of test utilities](#v8core-testutils), including functions for registering
tests with CMake; a tests runner, an extensive family of assertions and helper functions, and verbose error reporting.
Refer to the test toolkit's documentation for instructions on how to create and register tests.

Tests are placed in the `tests/v8` directory, grouped by modules.
Tests should also be topically grouped into test executables to reduce compilation times.
Each unit test function should test a single, or a small set of, API features, using the toolkit's assertion functions.
Particular attention should be given to edge cases and error conditions.
