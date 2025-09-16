# Geometry Functions

Functions to access geometric properties of the current cell and block
from pystencils kernels are available in the `sweepgen.symbolic` module.

## Blockforest, Block, and Cell Geometry

WaLBerla's blockforest comprises a hierarchy of geometry-defining objects:
The global simulation domain, which has one or more refinement levels,
each home to a set of blocks, which in turn define a grid of cells.
Many of these geometric properties are modelled as symbolic functions
and therefore accessible to kernels.
As Sweep kernels are always defined with respect to a single archetypical cell,
these symbolic functions always relate to that cell and its parent block.

The following tables lists the available symbolic functions,
the geometric concept they represent,
as well as equivalent methods of the waLBerla C++ API.
All Python APIs are listed relative to `sweepgen.symbolic`.

:::{admonition} Example Apps
:class: seealso

- {walberla-example-app}`DoubleShearLayer`
- {walberla-example-app}`ParallelPlates`
:::

### Domain Geometry

:::{list-table}
:header-rows: 1
:widths: auto

*   - Concept
    - Symbolic Functions
    - C++ API
*   - Domain Axis-Aligned Bounding Box
    - `domain.[x|y|z]_min()`, `domain.[x|y|z]_max()`
    - `StructuredBlockForest::getDomain()`
*   - Domain Cell Bounding Box (Current Refinement Level)
    - `domain_cell_bb.[x|y|z]_min()`, `domain_cell_bb.[x|y|z]_max()`
    - `StructuredBlockStorage::getDomainCellBB()`
:::

:::{note}
Similar to `CellInterval::max()` in the C++ API,
`domain_cell_bb.[x|y|z]_max()` refers to the largest valid cell index in the domain.
Therefore, a comparison with `<=` (instead of `<`) is necessary
in order to check if a cell is inside the domain.
:::

### Block Geometry

:::{list-table}
:header-rows: 1
:widths: auto

*   - Concept
    - Symbolic Functions
    - C++ API
*   - Block Axis-Aligned Bounding Box
    - `block.[x|y|z]_min()`, `block.[x|y|z]_max()`
    - `IBlock::getAABB()`
*   - Block Cell Bounding Box
    - `block_cell_bb.[x|y|z]_min()`, `block_cell_bb.[x|y|z]_max()`
    - `StructuredBlockForest::getBlockCellBB()`
:::

:::{note}
`block_cell_bb.[x|y|z]_max()` is *inclusive* to the block's cell bounding box;
to check if a cell lies inside the block, use `<=` as a comparison operator.
:::

### Cell Geometry

:::{list-table}
:header-rows: 1
:widths: auto

*   - Concept
    - Symbolic Functions
    - C++ API
*   - Cell Center Coordinates in Global Coordinate System
    - `cell.[x|y|z]()`
    - `StructuredBlockForest::getBlockLocalCellCenter()`
*   - Cell Center Coordinates in Block-Local Coordinate System
    - `cell.local_[x|y|z]()`
    - *No Equivalent*
*   - Cell Spacing
    - `cell.[dx|dy|dz]()`
    - `StructuredBlockForest::[dx|dy|dz]()`
*   - Cell Index in the Global Cell Grid
    - `cell_index.[x|y|z]_global()`
    - `StructuredBlockForest::transformBlockLocalToGlobalCell()`
*   - Cell Index in the Block-Local Cell Grid
    - `cell_index.[x|y|z]_local()`
    - `StructuredBlockForest::transformGlobalToBlockLocalCell()`
:::
