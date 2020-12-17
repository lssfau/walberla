************
Field module
************


Overview
========
The *field* contains classes and function exported from the waLBerla C++ field module.
*fields* are the central data structure in waLBerla where all simulation data is stored.
A *field* is basically a four dimensional array, where three dimensions are used to 
index cells in three dimensional space and the fourth dimension for storing multiple
values per cell. The field module supports two different layouts: array of structures (*zyxf*) and structure of
arrays (*fzyx*). 
For details have a look at the C++ documentation: :doxylink:class:`walberla::field::Field`



Fields and numpy
================

waLBerla *fields* can be easily converted to an  :py:class:`numpy.ndarray` using the array functionality of pybind11. This means
that a *numpy* array and a *field* can share the same memory such that fields can be manipulated using all the 
power of *numpy*::

   >>> import waLBerla
   >>> import numpy
   >>> field = waLBerla.field.createField( [3,3,3,1], float )
   >>> field[0,0,0,0]
   0.0
   >>> npArr = numpy.asarray( field )
   >>> npArr = waLBerla.field.toArray( field, True ) # convenience function, same as above (True includes ghostlayers)
   >>> npArr[:] = 42.0
   >>> field[0,0,0,0]
   42.0
   
A new *field* is created which is by default initialized with zero. Then a *numpy* array
is created which shares the same data. After modifying the *numpy* array also the field
values have changed.

A common source of error is to forget that some *numpy* functions create a copy of the data.
The copy is of course not shared with the field anymore::

   >>> npArr = numpy.roll( npArr, 1, axis=0 ) 
   >>> npArr[:] = 5
   >>> field[0,0,0]
   42.0

When during the array manipulation a copy was created the result has to be copied back into the
field again. Here the function :py:func:`numpy.copyto` is helpful:::
   
   >>> numpy.copyto( numpy.asarray( field ), npArr )
   >>> field = waLBerla.field.toArray( field ) # convenience function, equivalent to above
   >>> field[0,0,0]
   5.0 

   
Reference
=========


Classes
-------

.. py:class:: Field
   
   - Exported from C++ class :doxylink:class:`walberla::field::Field`
   - To modify or access a field class, the most convenient way is to create a *numpy.ndarray* view on it.
   
   .. py:method:: __array__
         
         The returned object implements pybind11::array and can be used for example to
         create a *numpy.ndarray* that shares the same data::
         
            numpy.asarray( field )

         With this function all ghostlayers are included in the view on the array.
         If the fourth dimension is one (this means only one value per cell). The returned numpy array has only 3 dimensions.
         
   
   .. py:method:: swapDataPointers ( otherField )
         
         Swaps the data of two fields. Only works if sizes, allocSizes and layout of the two fields
         are identical. The content of numpy arrays that have been created using the buffer interface
         are NOT swapped.
   

   The following attributes are read-only and provide information about sizes and memory layout 
   
   .. py:attribute:: size
   
         4-tuple with sizes of (x,y,z,f) coordinates not counting ghost layers 
         
   .. py:attribute:: allocSize
   
         The actual number of allocated elements for each coordinate. 
         Differences of size and allocSize are due to ghost layers and/or padding and depend on
         the chosen C++  :doxylink:class:`walberla::field::FieldAllocator` 
   
   .. py:attribute:: strides
         
         How many elements have to be skipped over in memory when incrementing the (x,y,z,f) dimension by one.  
   
   .. py:attribute:: offsets
      
         How many elements to skip over in memory from allocation begin to element (0,0,0,0)
   
   .. py:attribute:: layout
         
         Either *zyxf* (Array-of-Structures) or *fzyx* (Structure-of-Arrays) 
   
   
   
.. py:class:: GhostLayerField

   - Subclass of :py:class:`Field`
   - Exported from C++ class :doxylink:class:`walberla::field::GhostLayerField`

   .. py:attribute:: sizeWithGhostLayer
      
      4-tuple with sizes of (x,y,z,f) coordinates including ghost layers
   
   .. py:attribute:: nrOfGhostLayers
      
      The number of ghostlayers at each border of the field.

Free Functions
--------------

.. py:function:: createField( size, type, ghostLayers=1, layout = field.zyxf )

   Creates a new GhostLayerField
   
   :param size:        List of length 4 specifying x,y,z,f size of the field.
   :param type:        Type of the field elements. Valid types are the python types as well as some numpy types:
                        - Integer types: int, numpy.int[8,16,32,64]
                        - Unsigned types: numpy.uint[8,16,32,64]
                        - Float types : float, numpy.float32, numpy.float64
                        - Bool types : numpy.bool
                        
                       The type mapping is done via the C++ template trait ``walberla::python_coupling::isCppEqualToPythonType``
                       such that custom C++ types can be exported as well.
   :param ghostLayers: number of ghost layers of new field
   :param layout:      Either array-of-structures ``field.zyxf``  or structure-of-arrays  ``field.fzyx``



.. py:function:: addToStorage( blocks, name, dtype, fSize=1, layout=field.fzyx, ghostLayers=1, initValue=0.0, alignment=0)

   Adds a GhostLayerField to the given blockStorage
   
   :param blocks:       the structured blockstorage where the field should be added to
   :param name:         name of block data, is used to retrieve the created field later on
   :param dtype:        data type of the field
   :param fSize:        number of values per cell
   :param layout:       field.fzyx (SoA) or field.zyxf(AoS)
   :param ghostLayers:  number of ghost layers of the field
   :param initValue:    initial value for all cells, if None the types are default initialized (for most types zero)
   :param alignment:    alignment in bytes of the field vector


.. py:function:: gather( blocks, blockDataName, slice, targetRank=0 )
   
   Gathers part of the complete simulation domain (which is distributed to multiple processes)
   to one process. 
   
   :param blocks:        the blockstorage where the field is stored
   :param blockDataName: the name of the block data where the field is stored
   :param slice:         a slice object describing the region that should be collected in global coordinates
   :param targetRank:    world rank of process where the data should be gathered to
   
   Returns None on all processes except on process with rank targetRank. Here the collected field is returned.

   Slice gather example::
      >>> field.gather( blocks, 'density', makeSlice[ :,:,2 ] )
   
    
.. py:function:: createPackInfo( blocks, blockDataName )

   Creates a :doxylink:class:`walberla::field::communication::PackInfo` for a field.
   For details see tutorial on communication.
   
   
.. py:function:: createMPIDatatypeInfo( blocks, blockDataName )

   Creates a :doxylink:class:`walberla::field::communication::UniformMPIDatatypeInfo` for a field.
   For details see tutorial on communication.
                           