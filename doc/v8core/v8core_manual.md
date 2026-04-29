# V8 Core Library User Manual {#v8core-manual}

In the following pages, you find the user manual for the V8 core library.

## Using the V8 Core Library

### CMake

Activate the V8 core library in your build of waLBerla by setting the `WALBERLA_ENABLE_V8CORE` CMake option to `TRUE`;
e.g. on the command line:

```bash
cmake -DWALBERLA_ENABLE_V8CORE=TRUE
```

or in your CMake preset.
Then, in order to use the V8 core library in your application, link it against the `walberla::v8` library target:

```CMake
target_link_libraries( <YourApp> PRIVATE walberla::v8 )
```

### Include the Library Headers

In your C++ code, you can either bulk-include the entire V8 core library via `walberla/V8.hpp`,
or include its various submodules separately:

```
#include "walberla/V8.hpp"

// OR

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"
// ...
```

## Writing Portable Simulation Apps

The V8 core library is designed for *portability-first*, making it easy to write simulation apps
that can be compiled for, and run on, CPUs and GPU accelerators without changes to the code.
Portability in the V8 core is based on a twofold foundation:
 - Carefully designed library interfaces, using template metaprogramming protocols to distinguish
   between host and accelerator domains and hiding the details of the underlying parallel programming paradigm;
 - Automatic code generation using the SweepGen library for defining, parallelizing and offloading the heavy-duty numerical kernels.

Under the hood, waLBerla v8 uses plain C++20 plus OpenMP for parallelism on CPU,
and CUDA/HIP for targetting NVidia and AMD GPUs, respectively.

In the following, we explain
 - How to set up the source files and CMake targets for a portable simulation app;
 - How to use waLBerla's portability interfaces and the portable standard library.

### Sources and Build Setup

Create your application's primary translation unit (here `App.cpp`), and set up a CMake target:

```CMake
add_executable( App App.cpp )
```

When targetting CUDA or HIP, the app's source files must be interpreted accordingly as CUDA or HIP code.
To set the correct language, and pass the code to the correct compiler, call `walberla_set_gpu_language` on all
source files:

```CMake
walberla_set_gpu_language( App.cpp )
```

Finally, link against the `walberla::v8` library target:

```CMake
target_link_libraries( App PRIVATE walberla::v8 )
```

**Examples:**
 - @ref example-DoubleShearLayer

### Writing Device-Portable Code

Writing code (like numerical kernels) that should run on both CPU and GPU comes with a number of caveats.
Most importantly, such code must adhere to the restrictions imposed by the device programming language:
 - [Available C++ Features and Restrictions in CUDA device code](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/cpp-language-support.html)
 - [Available C++ Features and Restrictions in HIP device code](https://rocm.docs.amd.com/projects/HIP/en/latest/how-to/kernel_language_cpp_support.html)

CUDA and HIP device functions can only use APIs that are explicitly marked with the `__device__`-qualifier
(excepting `constexpr` functions, to a degree).

#### Device-Accessible waLBerla APIs

The V8 core library comes with a set of such device-accessible APIs.
Also, during its development, an increasing number of legacy APIs from the old `core` module are being made
available to the device.
You can find a list of waLBerla classes and functions available to device code at @ref v8core-device.

#### Portable Standard Library

For both CUDA and HIP there exist implementations of a select subset of the C++ standard library
ported to device code: [libcu++](https://nvidia.github.io/cccl/unstable/libcudacxx/index.html) (part of the CUDA SDK)
and [libhipcxx](https://github.com/ROCm/libhipcxx) (not yet part of the released HIP SDKs).
These libraries expose standard library classes, such as `std::array`, `std::span`, etc., to device-code.
The base namespaces for these two libraries are `cuda::std` and `hip::std`, respectively.

In waLBerla V8, we use these standard library ports in device-accessible code to facilitate portability;
for host-only code, we still rely on the default C++ standard library (`std::*`), and we recommend you do the same.

To transparently use the portable standard library, waLBerla defines the namespace alias `walberla::stdlib`,
which refers to either `cuda::std` or `hip::std` if CUDA or HIP are enabled,
or to the default `std` namespace if no GPU target is active.
To transparently include header files from the portable standard libraries,
use the `WALBERLA_STDLIB` macro.
Both `walberla::stdlib` and the `WALBERLA_STDLIB` include macro are available after including `walberla/V8.hpp`
or `walberla/v8/Device.hpp`:

```
#include "walberla/v8/Device.hpp"

//  Include headers from the portable standard library
#include WALBERLA_STDLIB(array)

//  Use portable classes:
walberla::stdlib::array< double, 3 >;
//  Refers to either `cuda::std::array` / `hip::std::array`,
//  or `std::array` if no GPU target is active
```
