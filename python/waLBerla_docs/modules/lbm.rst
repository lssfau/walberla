**********
LBM module
**********


Creation Functions
==================

.. py:function:: makeLatticeModel( stencil, collisionModel, forceModel, compressible, equilibriumAccuracyOrder=2 )

   Creates a new lattice model. A lattice model encapsulates all information about the lattice Boltzmann method.

   :param stencil:                  a string describing the stencil in DxQy notation e.g. 'D2Q9', 'D3Q19', 'D3Q27'
   :param collisionModel:           an instance of a collision model
   :param forceModel:               an instance of a force model
   :param compressible:             choose either a compressible or incompressible LBM scheme
   :param equilibriumAccuracyOrder: order of the equilibrium distribution. Valid values are 1 and 2. If not sure use 2 here.


   .. note ::
      The collision model and force model object are copied into the lattice model object. Changes to
      the passed force or collision model do not affect the state of the lattice model.
      Similarly after a sweep was created with a lattice model, the sweep itself is not changed when the lattice model
      it was created with was changed.



.. py:function:: addPdfFieldToStorage( blocks, name, latticeModel, initialVelocity=(0,0,0), initialDensity=1.0, ghostlayers=1, layout=field.zyxf, densityAdaptor="", velocityAdaptor="" )
                                       
   Adds a PDFField to the provided blockstorage and optionally a density and velocity adaptor.
   
   :param blocks:           blockstorage where the pdf field should be added to
   :param name:             block data id (string) of the new pdf field
   :param latticeModel:     see :py:meth:`makeLatticeModel` . The lattice model is copied into the pdf field. 
                            Later changes to the provided object do not affect the pdf field. To change parameters of 
                            the lattice model later, one has to iterate over all blocks, get the pdf field and retrieve  
                            a lattice model reference from it.
   :param initialVelocity:  lattice velocity the field is initialized with
   :param initialDensity:   density the field is initialized with
   :param ghostlayers:      number of ghost layers, has to be at least one
   :param layout:           memory layout of the field, ( see documentation of field module )
   :param densityAdaptor:   if a nonempty string is passed a :py:class:`FieldAdaptor` for the density is created with a 
                            blockdataID of the given name
   :param velocityAdaptor:  if a nonempty string is passed a :py:class:`FieldAdaptor` for the velocity is created with a 
                            blockdataID of the given name



      
.. py:function:: makeCellwiseSweep( blocks, pdfFieldID, flagFieldID="", flagList=[], velocityFieldID="" )

   Creates a new LBM sweep.
   
   :param blocks:          block storage where pdf field ( and if used, the flag field ) are stored
   :param pdfFieldID:      string identifier of the pdf field
   :param flagFieldID:     string identifier of the flag field. If empty string is passed, the LBM sweep is executed on all cells. 
   :param flagList:        Only necessary when a flagFieldID was specified. Pass a list of flag identifiers here,
                           describing the flags where the sweep should be executed
   :param velocityFieldID: optional velocity field ( field of fSize=3 and type=float) where the calculated velocity is written to.




.. py:class:: PdfField( field.GhostLayerField )
   
   .. py:attribute:: latticeModel:
   
   
   .. py:method:: setDensityAndVelocity( slice, velocity, density )
   
   .. py:method:: setToEquilibrium( slice, velocity, density )
   
   .. py:method:: getDensity( x,y,z )
   
   .. py:method:: getDensitySI ( x,y,z, rho_SI )   
   
   .. py:method:: getMomentumDensity ( x,y,z )   
   
   .. py:method:: getEquilibriumMomentumDensity ( x,y,z )   

   .. py:method:: getVelocity ( x,y,z )   
   
   .. py:method:: getVelocitySI ( x,y,z. dx_SI, dt_SI )
   
   .. py:method:: getEquilibriumVelocity ( x,y,z )   

   .. py:method:: getPressureTensor ( x,y,z )   
   
   


Boundary Handling
=================

.. py:class:: BoundaryHandling

   .. py:method:: isEmpty( x,y,z )
   
   .. py:method:: isNearBoundary( x,y,z )
   
   .. py:method:: isBoundary( x,y,z )
   
   .. py:method:: isDomain( x,y,z )
   
   .. py:method:: setDomain( x, y, z | slice )
   
   .. py:method:: forceDomain( x, y, z | slice )
   
   .. py:method:: fillWithDomain( x, y, z | slice | nrOfGhostLayersToInclude )
   
   .. py:method:: setBoundary( name, x, y, z | name, slice )
   
   .. py:method:: forceBoundary( name, x, y, z | name, slice )   
    
   .. py:method:: removeDomain( x, y, z | slice | nrOfGhostLayersToInclude )
      
   .. py:method:: removeBoundary( x, y, z | slice | nrOfGhostLayersToInclude )

   .. py:method:: clear( x, y, z | slice | nrOfGhostLayersToInclude )
      
      
Collision Models
================

.. py:class:: collisionModels.SRT

   Single Relaxation Time (BGK) lattice model   
   
   .. py:method:: __init__( omega, level=0 )

   
   .. py:attribute:: omega:
   
         Relaxation parameter ( = 1/tau )
   
   .. py:attribute:: viscosity:
   
   .. py:attribute:: level:
   
   .. py:method:: reset( omega, level=0 )
   
         Sets a new relaxation parameter for the given level


.. py:class:: collisionModels.SRTField( SRT )
   
   .. py:method:: __init__( omegaFieldID, level=0 )
   
      :param omegaFieldID:  this blockdata has to point to a floating point field of f-size=1 where for each cell
                            a different omega value is stored. 
                                  
   

.. py:class:: collisionModels.TRT
   
   .. py:method:: __init__( lambda_e, lambda_d, level=0 )
   
   .. staticmethod:: constructWithMagicNumber( omega, magicNumber=3.0/16.0 , level=0 )   
   
   .. py:attribute:: lambda_e:
   
   .. py:attribute:: lambda_d:
   
   .. py:attribute:: viscosity:

   .. py:attribute:: level:

   .. py:method:: reset( lambda_e, lambda_d, level=0 )
   
   .. py:method:: resetWithMagicNumber( omega, magicNumber=3.0/16.0 , level=0 )
   

.. py:class:: collisionModels.D3Q19MRT

   .. py:method:: __init__( s1, s2, s4, s9, s10, s16, level=0 )

   .. staticmethod:: constructTRTWithMagicNumber( omega, magicNumber=3.0/16.0 , level=0 )   
   
   .. staticmethod:: constructTRT( lambda_e, lambda_d, level=0 )
   
   .. staticmethod:: constructPanWithMagicNumber( omega, magicNumber=3.0/16.0 , level=0 )   
   
   .. staticmethod:: constructPan( lambda_e, lambda_d, level=0 )   
   
   .. py:attribute:: relaxationRates:
   
   .. py:attribute:: viscosity:

   .. py:attribute:: level:
   

      
Force Models
============

.. py:class:: forceModels.NoForce

.. py:class:: forceModels.SimpleConstant

   .. py:method:: __init__( force, level=0 )

.. py:class:: forceModels.EDMField

   .. py:method:: __init__( forceFieldID )

.. py:class:: forceModels.LuoConstant

   .. py:method:: __init__( force, level=0 )

.. py:class:: forceModels.LuoField

   .. py:method:: __init__( forceFieldID )

.. py:class:: forceModels.GuoConstant

   .. py:method:: __init__( force, level=0 )

.. py:class:: forceModels.Correction

   .. py:method:: __init__( previousMomentumDensityFieldID )
   
   
Lattice Models
==============

.. py:class:: LatticeModel

      Lattice models are created with the function :func:`makeLatticeModel` and encapsulate information about
      stencil, collision operator and force model.
      This information can be accessed through the following read-only attributes.

      .. py:attribute:: collisionModel:

            a *copy* of the collision model

      .. py:attribute:: forceModel:

            a *copy* of the force model

      .. py:attribute:: compressible:

            boolean signaling a compressible model

      .. py:attribute:: stencilName:

            a string describing the stencil in *DxQy* notation

      .. py:attribute:: communicationStencilName:

            name of stencil that should be used for communication. In most cases this is the same as stencilName

      .. py:attribute:: directions:

            a list of tuples containing the directions of the stencil. e.g. (0,0,0) for center, (1,0,0) for east etc.
            For a DxQy stencil the list as y entries, the tuples are of length x.




