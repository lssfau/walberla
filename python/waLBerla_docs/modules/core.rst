********************
waLBerla core Module
********************

Slices
======

Certain waLBerla functions require Python slices.
These slices can be created using the ``waLBerla.makeSlice`` object
with the usual Python list slicing syntax:

Example::

   import waLBerla
   waLBerla.field.gather(b, "SomeField" , waLBerla.makeSlice[:,0.5,:] )


Block Structure
===============

.. py:class:: CellInterval

   Export of :doxylink:class:`walberla::cell::CellInterval`.
   Describes a range of cells delimited by a min and a max cell. Caution: The max cell is included!

   .. py:staticmethod:: fromSlice( slice )

      Creates a CellInterval from a given Python slice e.g. ``CellInterval.fromSlice( makeSlice[0:2,0:1,0:4] )``

   .. py:method:: __init__( xMin, yMin, zMin, xMax, yMax, zMax )
      
      Constructs a cell interval using 6 integers corresponding to begin and end of the 
      CellInterval in each dimension. The maximum is included in the interval.
   
   .. py:attribute:: min
   
      Minimum cell.
   
   .. py:attribute:: max
   
      Maximum cell. 
   
   .. py:attribute:: size
   
      Size in x,y,z direction.
   
   .. py:attribute:: numCells
   
      Number of cells in the interval.
   
   .. py:method:: empty()
      
      True if CellInterval does not contain any cells.
   
   .. py:method:: positiveIndicesOnly()
   
      True if min and max cell are all non-negative.
   
   .. py:method:: contains( cell )
      
      Tests if a cell is contained in the CellInterval.
   
   .. py:method:: contains( cellInterval )
      
      Tests if another CellInterval is fully contained inside the CellInterval
   
   .. py:method:: overlaps( cellInterval )
   
      True if the CellIntervals have at least one common cell. 
   
   .. py:method:: shift( xShift, yShift, zShift )               
   
      Adds an offset to min and max cell, thus effectively shifting/moving the CellInterval.
   
   .. py:method:: getShifted( xShift, yShift, zShift )
   
      Returns a new CellInterval which is shifted by the given offsets
   
   .. py:method:: expand( nrOfCells )
   
      Subtracts nrOfCells from min, and adds nrOfCells to max. The size in on dimension
      effectively increases by 2*nrOfCells. 

   .. py:method:: getExpanded( x, y, z )

      Creates new CellInterval which is expanded by the given number of cells in each direction.

   .. py:method:: getExpanded( val )

      Shorthand for getExpanded( val, val, val )

   .. py:method:: intersect( cellInterval )
      
      Intersects the two CellIntervals and stores result in self.  
   
   .. py:method:: getIntersection( cellInterval )
      
      Returns the intersection interval between the current and the passed CellInterval. 
   

.. py:class:: AABB

   Export of :doxylink:class:`walberla::math::GenericAABB`.
   Axis aligned bounding box using floating point coordinates.
   
   .. py:method:: __init__( min, max )   
   
      Creates bounding box using a minimum and maximum point, given as tuples.
   
   .. py:method:: __init__( xMin, yMin, zMin, xMax, yMax, zMax )
   
      Creates bounding box using the given minimum and maximum values for each coordinate.
   
   .. py:attribute:: min
   
      Minimum cell.
   
   .. py:attribute:: max
   
      Maximum cell. 
   
   .. py:attribute:: size
   
      Size in x,y,z direction.
   
   .. py:method:: empty()
      
      True if enclosed volume is zero.
      
   .. py:method:: volume()
   
      Volume of the enclosed 3D cube.
      
   .. py:method:: center()
   
      Returns centroid of the enclosed cube.

      
   .. py:method:: contains( value )
      
      :param value: either another AABB or a point
      
   .. py:method:: containsClosedInterval( point, dx=0 )
   
      :param point: The point to be tested for containment
      :param dx:    An epsilon the box is extended by in each direction before the test
   
   .. py:method:: getExtended( x )
      
      :param x: either a scalar or 3-tuple with scaling factors for each dimension
      
   For the following methods see documentation of :doxylink:class:`walberla::math::GenericAABB`.
   
   .. py:method:: extend( scalarOrVector )
   
   .. py:method:: translate( vector )
   .. py:method:: getTranslated( translationVector )
   
   .. py:method:: scale( scalarOrVector )
   .. py:method:: getScaled( value )

   .. py:method:: merge( pointOrAABB )
   .. py:method:: getMerged( pointOrAABB )
   
   .. py:method:: intersect( aabb )
   .. py:method:: intersects( aabb, dx=0 )
   .. py:method:: intersectsClosedInterval( aabb, dx=0 )
   .. py:method:: intersectionVolume( aabb )
   .. py:method:: getIntersection( aabb )
   
   .. py:method:: isIdentical( aabb )
   .. py:method:: isEqual( aabb )
   
   .. py:method:: sqDistance( point )
   .. py:method:: sqSignedDistance( point )
   .. py:method:: sqMaxDistance( point )
   .. py:method:: distance( point )
   .. py:method:: signedDistance( point )
   .. py:method:: maxDistance( point )
   
  
  

.. py:class:: StructuredBlockForest
   
   StructuredBlockForest represents a collection of blocks. It can be created using the createUniformBlockGrid method.

   .. py:method:: getNumberOfLevels()
   .. py:method:: getDomain()
   
      Returns an axis aligned bounding box representing the complete simulation domain.
      
   .. py:method:: mapToPeriodicDomain( x,y,z )
   .. py:method:: mapToPeriodicDomain( point )
   .. py:method:: mapToPeriodicDomain( cell, level=0 )

   .. py:method:: getBlock( x,y,z )
   .. py:method:: containsGlobalBlockInformation( )
   .. py:method:: blocksOverlappedByAABB( point, aabb )
   .. py:method:: blocksContainedWithinAABB( point, aabb )
   
   .. py:method:: blockExists( point )
   .. py:method:: blockExistsLocally( point )
   .. py:method:: blockExistsRemotely( point )
   
   .. py:method:: atDomainXMinBorder( block )
   .. py:method:: atDomainYMinBorder( block )
   .. py:method:: atDomainZMinBorder( block )
   .. py:method:: atDomainXMaxBorder( block )
   .. py:method:: atDomainYMaxBorder( block )
   .. py:method:: atDomainZMaxBorder( block )

   .. py:method:: dx( level=0 )
   .. py:method:: dy( level=0 )
   .. py:method:: dz( level=0 )

   .. py:method:: getDomainCellBB( level=0 )
   .. py:method:: getBlockCellBB( block )
   .. py:method:: transformGlobalToLocal( block, object )
      
      :param object: either a cell (3 tuple) or a CellInterval
      
   .. py:method:: transformLocalToGlobal( block, object )

      :param object: either a cell (3 tuple) or a CellInterval
   
   
   .. py:attribute:: containsGlobalBlockInformation

   .. py:attribute:: periodic



.. py:class:: build_info

   .. py:attribute:: version

      Git Hash of waLBerla
      
   .. py:attribute:: type

      Type of build: Release,Debug, ...
      
   .. py:attribute:: compiler_flags
   .. py:attribute:: build_machine
   .. py:attribute:: source_dir
   .. py:attribute:: build_dir
         

Timing
======


.. py:class:: Timer

   .. py:method:: start()
   .. py:method:: stop()
   .. py:method:: reset()     
   .. py:method:: merge( otherTimer )
   
   .. py:attribute counter
   .. py:attribute total
   .. py:attribute sumOfSquares
   .. py:attribute average
   .. py:attribute variance
   .. py:attribute min      
   .. py:attribute max     
   .. py:attribute last      
      

.. py:class:: TimingPool
   
   .. py:method:: getReduced( reduceType, targetRank=0)
   
      :param reduceType: allowed values: total, min,avg, max 
       
   .. py:method:: merge( otherTimingPool, mergeDuplicates=True )
   
   .. py:method:: clear()
   .. py:method:: unifyRegisteredTimersAcrossProcesses()
   .. py:method:: logResultOnRoot()   
   

.. py:class:: TimingTree
   
   .. py:method:: start(timerName)

      :param value: name of the timer

   .. py:method:: stop(timerName)

      :param value: name of the timer

   .. py:method:: toDict()


Logging
=======

.. py:function:: abort( msg )


.. py:function:: log_devel( msg )


.. py:function:: log_devel_on_root( msg )


.. py:function:: log_result( msg )


.. py:function:: log_result_on_root( msg )


.. py:function:: log_warning( msg )


.. py:function:: log_warning_on_root( msg )


.. py:function:: log_info( msg )


.. py:function:: log_info_on_root( msg )


.. py:function:: log_progress( msg )


.. py:function:: log_progress_on_root( msg )


.. py:function:: log_detail( msg )


.. py:function:: log_detail_on_root( msg )




MPI
===

Rank & Communicator Infos
-------------------------


.. py:function:: rank()
.. py:function:: worldRank()
.. py:function:: numProcesses()
.. py:function:: hasCartesianSetup()
.. py:function:: rankValid( rank )

.. py:function:: worldBarrier()


Broadcast
---------

.. py:function:: broadcastInt( integerOrListOfIntegers, sendRank=0 )
.. py:function:: broadcastReal( realOrListOfReals, sendRank=0 )
.. py:function:: broadcastString( stringOrListOfStrings, sendRank=0 )



Reduction
---------

.. py:function:: reduceInt( integerOrListOfIntegers, operation, recvRank=0 )
.. py:function:: reduceReal( realOrListOfReals, operation, recvRank=0 )
.. py:function:: allreduceInt( integerOrListOfIntegers, operation )
.. py:function:: allreduceReal( realOrListOfReals, operation )


Gather
------

.. py:function:: gatherInt( integerOrListOfIntegers, recvRank=0 )
.. py:function:: gatherReal( realOrListOfReals, recvRank=0 )
.. py:function:: allgatherInt( integerOrListOfIntegers )
.. py:function:: allgatherReal( realOrListOfReals )


