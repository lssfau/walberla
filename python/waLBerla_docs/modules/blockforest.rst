******************
blockforest Module
******************



Reference
=========


.. py:function:: createUniformBlockGrid(blocks, cellsPerBlock, dx=1, oneBlockPerProcess=True, periodic=(0,0,0), keepGlobalBlockInformation=False)
   
   Creates a new uniform StructuredBlockForest. Similar to cpp function blockforest::createUniformBlockGrid.
   Specify blocks and cellsPerBlock.

   :param blocks:                       3-tuple with total number of blocks in x,y,z direction.
                                        When using this parameter you also have to pass cellsPerBlock.
   :param cellsPerBlock:                3-tuple with total number of cells per block in x,y,z direction.
                                        If this parameter is set, also blocks has to be set, but not cells
   :param dx:                           Side length of a single cell.
   :param oneBlockPerProcess:           If True, each process gets one block. If False, all blocks are put to one process.
                                        The second option makes only sense for debugging or testing.
   :param periodic:                     Periodicity of the domain in x,y,z direction
   :param keepGlobalBlockInformation:   If true, each process keeps information about remote blocks (blocks that reside
                                        on other processes). This information includes the process rank, the state, and
                                        the axis-aligned bounding box of any block (local or remote). [false by default]



   