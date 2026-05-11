# Boundary Handling

```{eval-rst}
.. module:: sweepgen.boundaries
```

SweepGen can generate two kinds of boundary handling sweeps for LBM simulations in waLBerla:
[sparse](#sparse_boundaries) and [grid-aligned](#grid_aligned_boundaries) boundary handlers.
Sparse boundary sweeps prescribe boundary conditions on arbitrary geometries by collecting
all lattice links crossing the boundary in an *index list*, and then iterating that list in parallel.
Grid-aligned boundary sweeps are more efficient, as they are implemented as a dense parallel iteration
over a cell interval, but can only handle grid-aligned boundary geometries with a fixed wall normal.
Their primary use case is boundary conditions at domain borders;
for most other cases, sparse boundaries should be used.

(sparse_boundaries)=
## Sparse Boundary Handlers

```{eval-rst}
.. autoclass:: sweepgen.boundaries.sparse.SparseBoundary
    :members:

.. autoclass:: sweepgen.boundaries.sparse.SparseBoundaryDefinition
    :members:

.. autoclass:: sweepgen.boundaries.sparse.FreeSlip
    :members:
```


(grid_aligned_boundaries)=
## Grid-Aligned Boundary Handlers

```{eval-rst}
.. autoclass:: sweepgen.boundaries.grid_aligned.GridAlignedNoSlip
    :members:

.. autoclass:: sweepgen.boundaries.grid_aligned.GridAlignedFreeSlip
    :members:
```
