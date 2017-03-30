******************
blockforest Module
******************



Reference
=========


.. py:function:: createUniformBlockGrid(cells, cellsPerBlock, blocks, periodic=(0,0,0), dx=1.0, oneBlockPerProcess=True)
   
   Creates a new uniform StructuredBlockStorage. Similar to cpp function createUniformBlockGridFromConfig.
   Specify either cells or (cellsPerBlock and blocks). 
   
   :param cells:              3-tuple with total numbers of cells in x,y,z direction. The returned BlockStorage may have
                              more cells if the cell count in a dimension is not divisible by the number of processes.
                              If this parameter is set, cellsPerBlock and blocks must not be set.
   :param cellsPerBlock:      3-tuple with total number of cells per block in x,y,z direction.
                              If this parameter is set, also blocks has to be set, but not cells
   :param blocks:             3-tuple with total number of blocks in x,y,z direction.
                              When using this parameter you also have to pass cellsPerBlock.
   :param periodic:           Periodicity of the domain in x,y,z direction
   :param dx:                 Side length of a single cell.
   :param oneBlockPerProcess: If True, each process gets one block. If False, all blocks are put to one process.
                              The second option makes only sense for debugging or testing.

   