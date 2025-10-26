# AMR Extraction Plan: Detailed Analysis

## Goal
Extract pure C++ adaptive mesh refinement (BlockForest) from waLBerla for FVM hydro development.

## Requirements
✅ BlockForest adaptive mesh system
✅ Field storage for FVM variables
✅ MPI support (optional - compile with/without)
✅ METIS support (optional - compile with/without)
✅ Structured grid (StructuredBlockForest)
✅ Load balancing (Hilbert curves + optional METIS)
✅ Dynamic refinement/coarsening
❌ No LBM, no particles, no GPU, no Python, no OpenMesh

---

## Module Structure to Extract

### PRIMARY MODULES (Must Extract)

#### 1. `/src/blockforest/` - **COMPLETE EXTRACTION** ✅
**Purpose**: Core adaptive mesh refinement engine

**Files (67 total):**

**Core AMR Files:**
- `BlockForest.h/.cpp` - Main AMR data structure
- `Block.h/.cpp` - Individual block representation
- `BlockID.h/.cpp` - Unique block identification
- `StructuredBlockForest.h/.cpp` - Structured grid wrapper
- `SetupBlockForest.h/.cpp` - Initial forest construction
- `SetupBlock.h/.cpp` - Setup phase block
- `PhantomBlockForest.h/.cpp` - Ghost forest for dynamic operations
- `PhantomBlock.h/.cpp` - Ghost blocks
- `BlockDataHandling.h/.cpp` - Data attachment to blocks
- `Types.h` - Type definitions
- `Utility.h` - Utility functions
- `all.h` - Convenience include

**Refinement Selection:**
- `AABBRefinementSelection.h` - Region-based refinement

**Initialization:**
- `Initialization.h/.cpp` - Forest initialization helpers
- `HilbertCurveConstruction.h` - Hilbert space-filling curves

**Reconstruction:**
- `BlockReconstruction.h/.cpp` - Block reconstruction from ID
- `BlockNeighborhoodConstruction.h/.cpp` - Neighbor relationships
- `BlockNeighborhoodSection.h` - Neighborhood sections

**Load Balancing:**
- `GlobalLoadBalancing.h` - Global load balancing algorithms

**Evaluation:**
- `BlockForestEvaluation.h/.cpp` - Statistics and analysis

**File I/O:**
- `BlockForestFile.h` - File format definitions
- `OutputColor.h/.cpp` - Visualization colors

**CMake:**
- `CMakeDefs.in.h` - CMake configuration header

**Subdirectories:**

`communication/` (9 files):
- `NonUniformBufferedScheme.h` - MPI communication for AMR
- `UniformBufferedScheme.h` - Uniform grid communication
- `UniformDirectScheme.h/.impl.h` - Direct communication
- `GeneratedNonUniformFieldPackInfo.h/.impl.h` - Field packing
- `NonUniformPackInfo.h` - Pack info interface
- `UniformToNonUniformPackInfoAdapter.h` - Adapter
- `DirectionBasedReduceScheme.h` - Reduction operations
- `LocalCommunicationMode.h` - Local communication modes
- `all.h`

`loadbalancing/` (14 files):
- `StaticCurve.h/.cpp` - Hilbert curve (NO METIS) ✅
- `DynamicCurve.h` - Dynamic Hilbert curve ✅
- `DynamicDiffusive.h` - Diffusive balancing ✅
- `Cartesian.h/.cpp` - Cartesian partitioning ✅
- `StaticParMetis.h/.cpp` - ParMETIS static (OPTIONAL - ifdef) ⚠️
- `DynamicParMetis.h/.cpp` - ParMETIS dynamic (OPTIONAL - ifdef) ⚠️
- `BlockInfo.h` - Block information
- `InfoCollection.h` - Info collection
- `NoPhantomData.h` - Phantom data interface
- `PODPhantomData.h` - POD phantom data
- `all.h`

`loadbalancing/level_determination/` (2 files):
- `MinMaxLevelDetermination.h/.cpp` - Min/max level constraints

`loadbalancing/weight_assignment/` (2 files):
- `WeightAssignmentFunctor.h` - Weight assignment ✅
- `MetisAssignmentFunctor.h` - METIS weights (OPTIONAL - ifdef) ⚠️

**Dependencies:**
- core/* (see below)
- domain_decomposition/*
- stencil/Directions.h, stencil/D3Q19.h

---

#### 2. `/src/field/` - **PARTIAL EXTRACTION**
**Purpose**: Field storage for FVM conservative variables (density, momentum, energy)

**Files to Extract (~20 files):**

**Core Field:**
- `Field.h/.impl.h` - Main field class ✅
- `GhostLayerField.h/.impl.h` - Ghost layer support ✅
- `AddToStorage.h` - Add fields to blocks ✅
- `FieldClone.h` - Field cloning ✅

**Iterators:**
- `iterators/FieldIterator.h` - Field iteration ✅
- `iterators/FieldPointer.h` - Field pointers ✅
- `iterators/IteratorMacros.h` - Iterator macros ✅

**Allocation:**
- `allocation/FieldAllocator.h` - Memory allocation ✅
- `allocation/FieldAlloc

atorStdVector.h` ✅

**Communication:**
- `communication/PackInfo.h` - MPI field packing ✅
- `communication/UniformMPIDatatypeInfo.h` - MPI datatypes ✅

**Refinement Support:**
- `refinement/PackInfo.h` - AMR field communication ✅
- `refinement/PackInfoHelper.h` - Refinement helpers ✅

**BlockForest Integration:**
- `blockforest/FieldAdaptor.h` - BlockForest integration ✅

**Optional but Useful:**
- `FlagField.h/.impl.h` - Boundary flag fields ✅
- `FileIO.h/.impl.h` - Field I/O ✅
- `Printers.h/.impl.h` - Debug printing ✅

**Skip:**
- vtk/ (we'll handle separately)
- python_coupling/
- cuda/
- opencl/
- Anything LBM-specific

**Dependencies:**
- core/*
- domain_decomposition/*
- stencil/*

---

#### 3. `/src/core/` - **PARTIAL EXTRACTION**
**Purpose**: Essential utilities and infrastructure

**Subdirectories to Extract:**

`math/` (COMPLETE):
- `AABB.h/.cpp` - Axis-aligned bounding boxes ✅
- `Vector3.h/.impl.h` - 3D vectors ✅
- `Matrix3.h` - 3x3 matrices ✅
- `Uint.h/.cpp` - Unsigned int utilities ✅
- `KahanSummation.h` - Accurate summation ✅
- `Random.h` - Random numbers ✅
- `IntegerFactorization.h` - Factorization ✅
- `Sample.h` - Sampling statistics ✅
- `DistributedSample.h` - MPI sampling ✅
- `FPClassify.h` - Floating point utilities ✅
- `Constants.h` - Math constants ✅

`mpi/` (COMPLETE if MPI enabled):
- `MPIManager.h/.cpp` - MPI initialization ✅
- `MPIWrapper.h` - MPI wrappers ✅
- `SendBuffer.h/.impl.h` - Send buffers ✅
- `RecvBuffer.h/.impl.h` - Receive buffers ✅
- `BufferSystem.h/.impl.h` - Buffer communication ✅
- `BufferSizeTrait.h` - Size traits ✅
- `Reduce.h/.impl.h` - MPI reductions ✅
- `Gatherv.h/.impl.h` - MPI gather ✅
- `SetReduction.h` - Set reductions ✅
- `Datatype.h` - MPI datatypes ✅
- `MPIHelper.h` - Helper functions ✅

`debug/` (COMPLETE):
- `Debug.h` - Debug macros ✅
- `CheckFunctions.h` - Assertion helpers ✅

`logging/` (COMPLETE):
- `Logging.h/.cpp` - Logging system ✅

`timing/` (COMPLETE):
- `Timer.h/.cpp` - Timing utilities ✅
- `TimingPool.h/.cpp` - Timing collections ✅

`uid/` (COMPLETE):
- `SUID.h/.cpp` - Set unique IDs ✅
- `UID.h` - Unique IDs ✅

`config/` (COMPLETE):
- `Config.h/.cpp` - Configuration file parser ✅

`cell/` (COMPLETE):
- `Cell.h` - Cell indices ✅
- `CellInterval.h` - Cell ranges ✅

`selectable/` (COMPLETE):
- `SetSelectableObject.h` - Selectable objects ✅
- `IsSetSelected.h` - Selection checking ✅

`load_balancing/` (PARTIAL - METIS optional):
- `MetisWrapper.h/.cpp` - METIS wrapper (OPTIONAL - ifdef) ⚠️
- `ParMetisWrapper.h/.cpp` - ParMETIS wrapper (OPTIONAL - ifdef) ⚠️

**Core Root Files:**
- `DataTypes.h` - Type definitions ✅
- `Abort.h/.cpp` - Abort functions ✅
- `Any.h` - std::any alternative ✅
- `Set.h` - Set utilities ✅
- `StringUtility.h` - String utils ✅
- `NonCopyable.h` - Non-copyable base ✅
- `OpenMP.h` - OpenMP helpers (OPTIONAL) ⚠️
- `EndianIndependentSerialization.h` - Serialization ✅
- `typeToString.h` - Type conversion ✅

**Skip:**
- python_coupling/
- cuda/
- waLBerlaBuildInfo.h (generated)
- Anything test-related

---

#### 4. `/src/domain_decomposition/` - **COMPLETE EXTRACTION**
**Purpose**: Block storage infrastructure

**Files (~15 files):**
- `IBlock.h/.cpp` - Block interface ✅
- `IBlockID.h` - Block ID interface ✅
- `BlockStorage.h/.cpp` - Block storage base ✅
- `StructuredBlockStorage.h/.cpp` - Structured storage ✅
- `BlockDataHandling.h` - Data handling interface ✅
- `MapPointToPeriodicDomain.h` - Periodic boundaries ✅
- `PeriodicIntersectionVolume.h` - Periodic intersections ✅
- `all.h` ✅

**Dependencies:**
- core/*
- stencil/*

---

#### 5. `/src/stencil/` - **PARTIAL EXTRACTION**
**Purpose**: Stencil definitions for FVM

**Files to Extract (~10 files):**
- `Directions.h` - Direction enumerations ✅
- `D3Q7.h` - 7-point stencil (Cartesian faces) ✅
- `D3Q19.h` - 19-point stencil (faces + edges) ✅
- `D3Q27.h` - 27-point stencil (full neighborhood) ✅
- `D2Q5.h, D2Q9.h` - 2D stencils ✅
- `Stencil.h` - Base stencil class ✅
- `Iterator.h` - Stencil iteration ✅

**Skip:**
- D3Q15, D3Q6, etc. (LBM-specific)

---

#### 6. `/src/boundary/` - **COMPLETE EXTRACTION** ✅ **REQUIRED FOR AMR**
**Purpose**: Boundary condition handling for AMR blocks (ESSENTIAL for FVM with ghost layers)

**Files to Extract (9 files - ALL):**
- `Boundary.h/.cpp` - Base boundary class ✅
- `BoundaryUID.h` - Boundary identifiers ✅
- `BoundaryHandling.h` - Main boundary manager (templated, 139KB) ✅
- `BoundaryHandlingCollection.h` - Multi-field boundaries ✅
- `ShiftedPeriodicity.h` - Periodic boundary support ✅
- `all.h` - Convenience include ✅
- `communication/HandlingPackInfo.h` - MPI boundary synchronization ✅
- `communication/all.h` - Convenience include ✅

**Why Required:**
- Per-block boundary flags (FlagField)
- Ghost layer communication for FVM stencils
- MPI synchronization of boundaries across processes
- Works naturally with AMR (each block has own boundaries)
- Foundation for implementing FVM boundary conditions

**Note**: This is generic infrastructure. You'll implement FVM-specific boundary conditions
(Dirichlet, Neumann, wall, inflow, outflow) by inheriting from `Boundary<>` template.

---

#### 7. `/src/timeloop/` - **PARTIAL EXTRACTION** (Optional but Useful)
**Purpose**: Time stepping infrastructure

**Files:**
- `TimeLoop.h/.cpp` - Main time loop ✅
- `ITimeStep.h` - Time step interface ✅
- `SweepTimeloop.h` - Sweep-based timeloop ✅

---

#### 8. `/src/vtk/` - **MINIMAL EXTRACTION** (Optional)
**Purpose**: VTK output for visualization

**Files:**
- `VTKOutput.h/.cpp` - VTK writer ✅
- `BlockCellDataWriter.h` - Cell data writer ✅

**Skip:**
- LBM-specific VTK writers
- Python VTK bindings

---

#### 9. `/src/geometry/` - **MINIMAL** (Optional)
**Purpose**: Basic geometric primitives

**Files:**
- Basic shape definitions if needed

---

### SUPPORTING FILES

#### CMake Build System

**Root:**
- `CMakeLists.txt` - Main CMake (HEAVILY MODIFIED) ⚠️
- `cmake/waLBerlaFunctions.cmake` - Helper functions
- `cmake/waLBerlaHelperFunctions.cmake` - More helpers
- `cmake/*.cmake` - Various utilities

**Each Module:**
- `src/blockforest/CMakeLists.txt` ✅
- `src/field/CMakeLists.txt` ✅
- `src/core/CMakeLists.txt` ✅
- etc.

**Configuration:**
- `waLBerlaDefinitions.h.in` - Main config header ⚠️
- `src/blockforest/CMakeDefs.in.h` ✅
- `src/field/CMakeDefs.in.h` ✅

---

## Dependency Graph

```
BlockForest
├── Block
│   ├── BlockID
│   ├── AABB (core/math)
│   ├── SUID (core/uid)
│   └── MPI buffers (core/mpi)
├── StructuredBlockForest
│   └── StructuredBlockStorage (domain_decomposition)
├── Communication
│   ├── BufferSystem (core/mpi)
│   └── Field PackInfo (field/communication)
├── LoadBalancing
│   ├── Hilbert curves (always)
│   ├── METIS (optional - ifdef)
│   └── ParMETIS (optional - ifdef)
└── Domain Decomposition
    ├── BlockStorage
    ├── IBlock
    └── Stencil (stencil/)

Field
├── Field<T> template
├── GhostLayerField
├── BlockForest integration (field/blockforest)
└── Refinement (field/refinement)

Core
├── DataTypes
├── Math (AABB, Vector3, Uint)
├── MPI (optional - compile with/without)
├── Logging
├── Timing
├── Config
└── METIS wrappers (optional - ifdef)
```

---

## Compilation Modes

### Mode 1: Serial (No MPI, No METIS)
```cmake
-DWALBERLA_BUILD_WITH_MPI=OFF
-DWALBERLA_BUILD_WITH_METIS=OFF
```
- Hilbert curve load balancing only
- Single process
- Minimal dependencies

### Mode 2: MPI (No METIS)
```cmake
-DWALBERLA_BUILD_WITH_MPI=ON
-DWALBERLA_BUILD_WITH_METIS=OFF
```
- Hilbert curve + Diffusive load balancing
- Multi-process
- Requires: MPI library

### Mode 3: MPI + METIS
```cmake
-DWALBERLA_BUILD_WITH_MPI=ON
-DWALBERLA_BUILD_WITH_METIS=ON
-DWALBERLA_BUILD_WITH_PARMETIS=ON
```
- All load balancing options
- Multi-process
- Optimal partitioning
- Requires: MPI, METIS, ParMETIS

---

## External Dependencies

### Required:
- **CMake >= 3.26**
- **C++20 compiler** (GCC 10+, Clang 10+, MSVC 2019+)

### Optional:
- **MPI library** (OpenMPI, MPICH, Intel MPI) - for parallel
- **METIS** - for optimal load balancing
- **ParMETIS** - for parallel optimal load balancing

---

## File Count Estimate

| Module | Files | Lines (est.) |
|--------|-------|--------------|
| blockforest/ | 67 | 15,000 |
| field/ | 20 | 8,000 |
| core/ | 80 | 25,000 |
| domain_decomposition/ | 15 | 5,000 |
| stencil/ | 10 | 2,000 |
| boundary/ | 9 | 5,000 |
| communication/ | 5 | 2,000 |
| timeloop/ | 5 | 1,000 |
| vtk/ | 3 | 1,000 |
| CMake files | 20 | 2,000 |
| **TOTAL** | **~234** | **~66,000** |

---

## Extraction Strategy

### Phase 1: Create New Directory Structure
```
walberla_amr/
├── src/
│   ├── blockforest/      ✅ AMR core
│   ├── field/            ✅ Field storage
│   ├── core/             ✅ Utilities
│   ├── domain_decomposition/ ✅ Block infrastructure
│   ├── stencil/          ✅ Stencil definitions
│   ├── boundary/         ✅ Boundary handling (REQUIRED)
│   ├── communication/    ✅ MPI communication
│   ├── timeloop/         ✅ Time stepping
│   └── vtk/              ✅ Visualization
├── cmake/
├── CMakeLists.txt
├── examples/
│   └── simple_amr_fvm_test.cpp
└── README.md
```

### Phase 2: Copy Files Module-by-Module
1. Copy blockforest/ complete
2. Copy required core/ subdirectories
3. Copy domain_decomposition/ complete
4. Copy field/ (filtered)
5. Copy stencil/ (filtered)
6. Copy supporting modules

### Phase 3: CMake Configuration
1. Create minimal root CMakeLists.txt
2. Optional MPI detection
3. Optional METIS detection
4. Module-specific CMakeLists.txt
5. Generate configuration headers

### Phase 4: Testing
1. Test compilation: Serial mode
2. Test compilation: MPI mode
3. Test compilation: MPI+METIS mode
4. Create simple AMR test

---

## Questions Before Proceeding

1. **Directory structure**: Should I extract to a new subdirectory in the current repo, or create a completely separate clean directory?

2. **VTK output**: Do you want VTK file output for visualization, or will you handle that separately?

3. **Timeloop**: Do you want the time-stepping infrastructure, or will you write your own time loop for FVM?

4. **Boundary conditions**: Do you want waLBerla's boundary handling, or will you implement FVM-specific boundaries yourself?

5. **Examples**: Should I create a simple AMR test example as part of the extraction?

---

## Risks & Mitigation

| Risk | Mitigation |
|------|------------|
| Missing hidden dependencies | Systematic grep through all includes |
| METIS ifdef guards broken | Test all 3 compilation modes |
| MPI ifdef guards broken | Test with/without MPI |
| CMake configuration issues | Minimal CMake, test incrementally |
| Template instantiation errors | Include .impl.h files properly |

---

## Estimated Time

- **Planning & Approval**: 2 hours ✅
- **File extraction**: 3-4 hours
- **CMake setup**: 2-3 hours
- **Testing & debugging**: 3-5 hours
- **Example creation**: 1-2 hours
- **TOTAL**: 11-16 hours

---

## FINAL EXTRACTION SET (Updated with Boundary Module)

**Complete extraction includes:**
```
✅ blockforest/ - Complete (67 files) - AMR core engine
✅ field/ - Field storage (20 files) - for FVM conservative variables
✅ core/ - Essential utilities (80 files) - math, MPI, timing, logging, config
✅ domain_decomposition/ - Complete (15 files) - block storage infrastructure
✅ stencil/ - Stencil definitions (10 files) - D3Q7, D3Q19, D3Q27 for FVM
✅ boundary/ - Complete (9 files) - REQUIRED for AMR boundary handling
✅ communication/ - MPI wrappers (5 files) - parallel communication
✅ timeloop/ - Time stepping (5 files) - useful for FVM time integration
✅ vtk/ - Visualization (3 files) - VTK output for debugging

Total: ~234 files, ~66,000 lines
```

**Compilation modes:**
- Serial (no MPI, no METIS) - single process, Hilbert curves
- MPI (no METIS) - parallel, Hilbert + Diffusive balancing
- MPI + METIS (optional) - parallel, optimal graph partitioning

**CMake defaults:**
- MPI: ON (compile with/without via flag)
- METIS: OFF (enable optionally)
- OpenMP: OFF (not needed initially)

---

## Ready to Proceed?

Please review this plan and confirm:
1. **Directory structure**: Extract to new `walberla_amr/` subdirectory? Or separate location?
2. **Module selection**: All modules above OK? (boundary now included ✅)
3. **Any concerns or modifications?**

Once approved, I will execute systematically and carefully.
