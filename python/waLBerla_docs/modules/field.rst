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

waLBerla *fields* can be easily converted to an  :py:class:`numpy.ndarray` using Python's buffer protocol. This means
that a *numpy* array and a *field* can share the same memory such that fields can be manipulated using all the 
power of *numpy*::

   >>> import waLBerla
   >>> import numpy
   >>> field = waLBerla.field.createField( [3,3,3,1], float )
   >>> field[0,0,0,0]
   0.0
   >>> npArr = numpy.asarray( field.buffer() )
   >>> npArr = waLBerla.field.toArray( field ) # convenience function, same as above
   >>> npArr[:] = 42.0
   >>> field[0,0,0,0]
   42.0
   
A new *field* is created which is by default initialized with zero. Then a *numpy* array
is created which shares the same data. After modifying the *numpy* array also the field
values have changed. To view the field including ghost layers additional parameters to 
:py:meth:`Field.buffer` are required.

A common source of error is to forget that some *numpy* functions create a copy of the data.
The copy is of course not shared with the field anymore::

   >>> npArr = numpy.roll( npArr, 1, axis=0 ) 
   >>> npArr[:] = 5
   >>> field[0,0,0]
   42.0

When during the array manipulation a copy was created the result has to be copied back into the
field again. Here the function :py:func:`numpy.copyto` is helpful:::
   
   >>> numpy.copyto( numpy.asarray( field.buffer() ), npArr )
   >>> field = waLBerla.field.fromArray( npArr ) # convenience function, equivalent to above
   >>> field[0,0,0]
   5.0 

   
Reference
=========


Classes
-------

.. py:class:: Field
   
   - Exported from C++ class :doxylink:class:`walberla::field::Field`
   - To modify or access a field class, the most convenient way is to create a *numpy.ndarray* view on it.
   
   .. py:method:: buffer ( withGhostLayers=False )
         
         The returned object implements the Python Buffer Protocol and can be used for example to
         create a *numpy.ndarray* that shares the same data::
         
            numpy.asarray( field.buffer(withGhostLayers=True) )
         
         The optional parameter specifies if the ghost layers are part of the buffer.
         
   
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
   

.. py:class:: FlagField
   
   Subclass of :py:class:`GhostLayerField` where the value type is an unsigned integer and the
   size of the f coordinate is fixed to one element. FlagFields provide additional management function
   for storing multiple booleans per cell (encoded in bits). 
   
   
   FlagFields are exported from C++ class :doxylink:class:`walberla::field::FlagField`
   
   .. py:method:: registerFlag( flagName, bitNr = None )
   
         Reserves the next free bit (if bitNr is None ) or the specified bit using the provided flag name.
         Returns an integer where the reserved bit is set to one, all other bits are set to zero.
                  
   .. py:method:: flag( flagname )
      
         Returns an integer where the specified flag is set to one, all other flags are zero.
      
   .. py:method:: flagName( flag )
         
         Maps from integer where on bit is set to the name of the flag.
      
   .. py:attribute:: flags
          
          List with registered flag names.
          


.. py:class:: FieldAdaptor

   A field adaptor is an object that emulates a GhostLayerField but does not store data itself.
   Adaptors can only be created by C++ using :doxylink:class:`walberla::field::GhostLayerFieldAdaptor`.
   
   When accessing a cell of an adaptor, its value is computed on the fly based on one or multiple input fields.
   A VelocityAdaptor, for example, computes the macroscopic velocity in a cell based on a field of particle distribution functions (PDFs).
   Since adaptor do not hold data themselves they cannot be converted directly to numpy arrays ( see :py:meth:`copyToField` ).
   Since this operation is expensive consider accessing only the required adaptor values using the getitem operator. 
   
   .. py:method:: copyToField ()
   
      Creates a field by computing the adaptor value for every cell (potentially expensive). 
      Returns this temporary field. Modifications of this field
      do not affect the adaptor or the adaptor base field. 
   
    
   .. py:attribute:: size
   
         4-tuple with sizes of (x,y,z,f) coordinates not counting ghost layers 
    
   .. py:attribute:: sizeWithGhostLayer
      
      4-tuple with sizes of (x,y,z,f) coordinates including ghost layers
    



Free Functions
--------------

.. py:function:: createField( size, type, ghostLayers=1, layout = field.zyxf )

   Creates a new GhostLayerField
   
   :param size:        List of length 3 or 4 specifying x,y,z,f size of the field. 
                       If list is of length 3 f-size is assumed to be 1
   :param type:        Type of the field elements. Valid types are the python types as well as some numpy types:
                        - Integer types: int, numpy.int[8,16,32,64]
                        - Unsigned types: numpy.uint[8,16,32,64]
                        - Float types : float, numpy.float32, numpy.float64
                        
                       The type mapping is done via the C++ template trait ``walberla::python_coupling::isCppEqualToPythonType``
                       such that custom C++ types can be exported as well.
   :param ghostLayers: number of ghost layers of new field
   :param layout:       Either array-of-structures ``field.zyxf``  or structure-of-arrays  ``field.fzyx``
   
              
       
             

.. py:function:: createFlagField( size, nrOfBits=32, ghostLayers=1 )
   
   Creates a new FlagField
   
   :param size:        list of length 3 with x,y,z size of field
   :param nrOfBits:    how many flags can be stored per cell. Allowed values are 8,16,32,64
   :param ghostLayers: number of ghost layers of new field



.. note:: The ValueError "Cannot create field of this (type,f-size) combination"
          means that in C++ this specific choice of type and f-size was not exported to Python.
          In C++ these are template parameters, so a separate field class has to be instantiated for each
          combination. 




.. py:function:: addToStorage( blocks, name, type, fSize=1, ghostLayers=1, layout=field.fzyx, initValue=None)

   Adds a GhostLayerField to the given blockStorage
   
   :param blocks:    the structured blockstorage where the field should be added to
   :param name:      name of block data, is used to retrieve the created field later on
   :param initValue: initial value for all cells, if None the types are default initialized (for most types zero)
   
   The remaining parameter are the same as in  :py:func:`createField`
    
   
   
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
                           