# Timing Pool {#example-TimingPool}

\viewOnGitlab{apps/examples/TimingPool}

\tableofcontents

## Domain Setup and Initial State

To showcase the use of waLBerla's Timing Pools,
a trivial 3D domain is setup as a arbitrary cuboid with uniform spacing in every direction.
The workload distribution is determined by waLBerla according to the available MPI-processes.

This example is not actually solving a physical problem but simply passing some time in bogus sweeps to measure some timings.
It defines three Sweeps with different idling times as well as some synchronization between them to showcase the runtime distribution between Sweeps.

WaLBerla supports three different kinds of Timers:

- `WcTimer`: **Wall-clock Timer.** Measures the actual elapsed real-world time. This includes time spent waiting for I/O, system delays, or other processes. It is based on `std::chrono::high_resolution_clock`.
- `CpuTimer`: **CPU Execution Timer.** Measures the amount of time the processor spent actively executing the process's instructions in user mode (using `getrusage`). This excludes time when the process is suspended or waiting for external resources.
- `DeviceSynchronizeTimer`: **Blocking GPU-Aware Wall-clock Timer.** Operates like the `WcTimer` but explicitly invokes a device synchronization (e.g., `gpuDeviceSynchronize`) before stopping. This ensures that all asynchronous GPU kernels have finished execution. But comes to the cost of parallel host- executions while waiting for the work on the device.

Every Timing is given in seconds \[s\].

## Build and config

To run the Timing Pool example application waLBerla must be build with `WALBERLA_BUILD_EXAMPLES=true`.
There are no other general setup restrictions for this example.

## Application Frame

The C++ Example application is implemented in `TimingPoolExample.cpp`.

### Preamble

To inspect Sweep based timer statistics from the Timing Pool,
the TimingPool header must be included:

```cpp
#include "core/timing/TimingPool.h"
```

### General Scope

The execution structure within the `run` method is common for most waLBerla applications.

1. Initialization of the simulations's environment.
2. Extracting a `Config` object from the environment.
3. Initializing the `BlockForest`.
4. Registering the simulations `Field`s onto the `BlockForest`.
5. Setting up the communication pattern.
6. Extracting the simulations runtime parameters.
7. **Defining the simulation's Execution Graph, i.e., registering `Sweep`s to the `TimeLoop`.**
8. Preparing the runtime watchers/observers/metrics/measure-tools, e.g., the timers, and the `TimingPool`.
9. Performing the simulation run.
10. Finalizing the Simulation.

- A `TimingPools` is an optional parameter for executing a pre-registered `TimingLoop`.
- Every Sweep is an element in the `TimingPool`.

> [!note]
> This documentation only discusses codesections relevant to understand the use of TimingPools.
> For guidance regarding the simulation setup,
> please refer to [Basic Fluid Simulations](#example-basic-fluid-simulations)
> or the [Tutorials](#tutorials) section.

### Parameters and Initial State

The runtime distribution of each Sweeps can be controlled by adjusting the `delayFactorKernel` parameters specified in [TimingPoolExample.prm](./TimingPoolExample.prm):

\dontinclude[strip] TimingPool/TimingPoolExample.prm
\skip Parameters
\until }

> [!note]
> The parameter `gatherStatistics` switches to a more complicated runtime distribution and can be ignored for now.

### Simulation Loop

The execution flow of waLBerla applications is defined in `TimingLoop`s where different executions Sweeps are attached to.

\dontinclude[strip] TimingPool/TimingPoolExample.cpp
\skip Begin-Timing-Pool-Definition
\until End-Timing-Pool-Definition

### Prepare Time Measurement

Time Measurements in waLBera make use of two timing objects.

#### `WCTimer`

The `WCTimer` object works as simple stopwatch for program sections,
by providing a wrapper for the
[C++ chrono library](https://en.cppreference.com/w/cpp/chrono.html#chrono_library)
and waLBerla's MPI Module.

#### `WCTimingPools`

The `WCTimingPool` object is a collection of Sweep-specific `WCTimer`s.
When registered for a `TimingLoop` execution, they gather timing information of every individual execution Sweep within the Loop.
This way it is possible to investigate the runtime distribution of a simulation to find bottle necks.

The TimingPools creates min, max and average statistics per timestep and thread, facilitating the computation of performance metrics and runtime analysis.

### Invoke the Timed Timeloop

\snippet TimingPool/TimingPoolExample.cpp Marker-Timing-Pool-Execution

### Evaluate Runtime measurements

The `TimingPool` provides runtime reduction over all timestops.
The thread specific runtime reduction is handled by the user.

Since our LBM based simulations usually require synchronization between every timestep,
we are usually interested in the actually time we had to wait for the result.
So the time it took for the slowest thread to execute the kernel.
This can be obtained by a reducing with the `Max` operator.

```cpp
double time = simTimer.max();
WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
```
