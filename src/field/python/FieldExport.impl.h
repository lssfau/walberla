//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file FieldExport.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"


#include "core/logging/Logging.h"
#include "core/VectorTrait.h"
#include "field/Field.h"
#include "field/GhostLayerField.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"

#include "field/AddToStorage.h"
#include "field/python/GatherExport.h"
#include "field/vtk/VTKWriter.h"
#include "field/vtk/FlagFieldMapping.h"

#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"

#include <boost/mpl/vector.hpp>

#include <iostream>
#include <type_traits>

namespace walberla {
namespace field {



namespace internal {

   //===================================================================================================================
   //
   //  Buffer Protocol for Fields
   //
   //===================================================================================================================



   /* This section implements the Python buffer protocol for walberla::Field's
    *
    *  - Why? The buffer protocol enables other Python types to use the memory belonging to a field.
    *         One can for example construct a numpy.array that operates on the field data
    *  - How? a = numpy.asarray( myWalberlaField.buffer() )
    *         creates a numpy array which uses the field data in a read-write way! no data is copied
    *  - Python buffer protocol: http://docs.python.org/dev/c-api/buffer.html
    *  - Why so complicated?
    *       boost::python does not yet (as in version 1.55) support the buffer protocol
    *       so everything has to be written in the native Python C interface.
    *       In order to export the Field with boost python and keep the native C interface part as
    *       small as possible, a new Type called FieldBuffer is introduced which is exported in the
    *       C interface. This type can only be created using the buffer() function of the field.
    *
    *
    *  - Lifetime issues:
    *       0) f   = walberla.create_field( ( 5,3,2) )
    *       1) buf = f.buffer():
    *             - creates a FieldBuffer object which has as only member the Field
    *             - f is not deallocated as long as 'buf' exists
    *       2) a = numpy.asarray( buf )
    *             - calls 'fieldbuffer_get' which creates a Py_buffer,
    *                extracts the information, then immediately calls 'fieldbuffer_release'
    *                to clean up the Py_buffer.
    *             - buf still exists (   fieldbuffer_release only cleans up Py_buffer not
    *                the FieldBuffer object 'buf' )
    *       3) del f
    *          del b
    *             - the buffer object b is not deallocated since a still has a reference to it
    *             - since b is not deleted f is not deleted since b has still a reference to it
    *               i.e   a -> b -> f
    *       4) when a is deleted everything can be cleaned up, which means that fieldbuffer_dealloc
    *          is called, which decrements the reference count for f, so that f can be deallocated if
    *          not used otherwise
    *
    */
   template<class T> struct PythonFormatString                    { inline static char * get() { static char value [] = "B"; return value; } };

   template<>        struct PythonFormatString<double>            { inline static char * get() { static char value [] = "d"; return value; } };
   template<>        struct PythonFormatString<float>             { inline static char * get() { static char value [] = "f"; return value; } };
   template<>        struct PythonFormatString<unsigned short>    { inline static char * get() { static char value [] = "H"; return value; } };
   template<>        struct PythonFormatString<int>               { inline static char * get() { static char value [] = "i"; return value; } };
   template<>        struct PythonFormatString<unsigned int>      { inline static char * get() { static char value [] = "I"; return value; } };
   template<>        struct PythonFormatString<long>              { inline static char * get() { static char value [] = "l"; return value; } };
   template<>        struct PythonFormatString<unsigned long>     { inline static char * get() { static char value [] = "L"; return value; } };
   template<>        struct PythonFormatString<long long>         { inline static char * get() { static char value [] = "q"; return value; } };
   template<>        struct PythonFormatString<unsigned long long>{ inline static char * get() { static char value [] = "Q"; return value; } };


   typedef struct {
       PyObject_HEAD
       PyObject * field;
       PyObject * fieldAlloc;
       void * fieldData;
   } FieldBufferPyObject;


   template<typename T, uint_t fs>
   struct FieldBufferGetDispatch
   {
        static int get( FieldBufferPyObject * exporter, Py_buffer * view, int flags, bool withGhostLayer )
        {
           namespace bp = boost::python;

           bp::object fieldObject ( bp::handle<>(bp::borrowed( exporter->field ) ) );
           Field<T,fs> * field = bp::extract< Field<T,fs> * > ( fieldObject );

           bp::object fieldAllocObject ( bp::handle<>(bp::borrowed( exporter->fieldAlloc ) ) );
           FieldAllocator<T> * fieldAlloc = bp::extract< FieldAllocator<T> * > ( fieldAllocObject );
           fieldAlloc->incrementReferenceCount( field->data() );

           view->obj = (PyObject*) exporter;
           Py_INCREF( view->obj );

           uint_t size[3];
           auto glField = dynamic_cast<GhostLayerField<T,fs> * >( field );
           if ( glField && withGhostLayer )
           {
              size[0] = glField->xSizeWithGhostLayer();
              size[1] = glField->ySizeWithGhostLayer();
              size[2] = glField->zSizeWithGhostLayer();
              cell_idx_t gl = cell_idx_c( glField->nrOfGhostLayers() );
              view->buf = &( field->get( -gl, -gl, -gl,0 ) );
           }
           else
           {
              size[0] = field->xSize();
              size[1] = field->ySize();
              size[2] = field->zSize();
              view->buf = & ( field->get(0,0,0,0) ) ;
           }


           // Mandatory
           view->len = Py_ssize_t( size[0] * size[1] * size[2] * fs * sizeof(T) );
           view->itemsize = sizeof( T );

           view->ndim = 4;

           view->format = NULL;
           if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
              view->format = PythonFormatString<T>::get();

           view->shape=NULL;
           if ((flags & PyBUF_ND) == PyBUF_ND)
           {
              view->shape = new Py_ssize_t[4];
              view->shape[0u] = Py_ssize_t( size[0u] );
              view->shape[1u] = Py_ssize_t( size[1u] );
              view->shape[2u] = Py_ssize_t( size[2u] );
              view->shape[3u] = Py_ssize_t( field->fSize() );
           }

           view->strides = NULL;
           if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
           {
              view->strides = new Py_ssize_t[4u];
              view->strides[0u] = Py_ssize_t( uint_c( field->xStride() ) * sizeof(T) );
              view->strides[1u] = Py_ssize_t( uint_c( field->yStride() ) * sizeof(T) );
              view->strides[2u] = Py_ssize_t( uint_c( field->zStride() ) * sizeof(T) );
              view->strides[3u] = Py_ssize_t( uint_c( field->fStride() ) * sizeof(T) );
           }

           view->suboffsets = NULL;
           view->internal = NULL;
           view->readonly = false;

           // We don't need the field any more. For freeing only the allocator and the field data is necessary
           Py_DECREF( exporter->field );
           exporter->field = NULL;

           return 0;
        }
    };

    template<typename VectorType>
    struct FieldBufferGetDispatch<VectorType,1>
    {
        static int get( FieldBufferPyObject * exporter, Py_buffer * view, int flags, bool withGhostLayer )
        {
           namespace bp = boost::python;

           typedef VectorTrait<VectorType> VecTrait;
           typedef typename VecTrait::OutputType ElementType;

           static_assert( sizeof(VectorType) == VecTrait::F_SIZE*sizeof(ElementType),
                          "Creating Python Memory View works only for vectors types that hold their elements consecutively in memory" );

           bp::object fieldObject ( bp::handle<>(bp::borrowed( exporter->field ) ) );
           Field<VectorType,1> * field = bp::extract< Field<VectorType,1> * > ( fieldObject );

           bp::object fieldAllocObject ( bp::handle<>(bp::borrowed( exporter->fieldAlloc ) ) );
           FieldAllocator<VectorType> * fieldAlloc = bp::extract< FieldAllocator<VectorType> * > ( fieldAllocObject );
           fieldAlloc->incrementReferenceCount( field->data() );

           view->obj = (PyObject*) exporter;
           Py_INCREF( view->obj );

           uint_t size[3];
           auto glField = dynamic_cast<GhostLayerField<VectorType,1> * >( field );
           if ( glField && withGhostLayer )
           {
              size[0] = glField->xSizeWithGhostLayer();
              size[1] = glField->ySizeWithGhostLayer();
              size[2] = glField->zSizeWithGhostLayer();
              cell_idx_t gl = cell_idx_c( glField->nrOfGhostLayers() );
              view->buf = &( field->get( -gl, -gl, -gl,0 ) );
           }
           else
           {
              size[0] = field->xSize();
              size[1] = field->ySize();
              size[2] = field->zSize();
              view->buf = & ( field->get(0,0,0,0) ) ;
           }


           // Mandatory
           view->len = Py_ssize_t( size[0] * size[1] * size[2] * VecTrait::F_SIZE * sizeof(ElementType) );
           view->itemsize = sizeof( ElementType );

           view->ndim = 4;

           view->format = NULL;
           if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
              view->format = PythonFormatString<ElementType>::get();

           view->shape=NULL;
           if ((flags & PyBUF_ND) == PyBUF_ND)
           {
              view->shape = new Py_ssize_t[4];
              view->shape[0u] = Py_ssize_t( size[0u] );
              view->shape[1u] = Py_ssize_t( size[1u] );
              view->shape[2u] = Py_ssize_t( size[2u] );
              view->shape[3u] = Py_ssize_t( VecTrait::F_SIZE );
           }

           view->strides = NULL;
           if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
           {
              view->strides = new Py_ssize_t[4u];
              view->strides[0u] = Py_ssize_t( uint_c( field->xStride() ) * VecTrait::F_SIZE * sizeof(ElementType) );
              view->strides[1u] = Py_ssize_t( uint_c( field->yStride() ) * VecTrait::F_SIZE * sizeof(ElementType) );
              view->strides[2u] = Py_ssize_t( uint_c( field->zStride() ) * VecTrait::F_SIZE * sizeof(ElementType) );
              view->strides[3u] = Py_ssize_t( sizeof(ElementType) );
           }

           view->suboffsets = NULL;
           view->internal = NULL;
           view->readonly = false;

           // We don't need the field any more. For freeing only the allocator and the field data is necessary
           Py_DECREF( exporter->field );
           exporter->field = NULL;

           return 0;
        }
    };

   template<typename T, uint_t fs>
   int fieldbuffer_get_withGl ( FieldBufferPyObject * exporter, Py_buffer * view, int flags )
   {
      return FieldBufferGetDispatch<T,fs>::get( exporter, view, flags, true );
   }

   template<typename T, uint_t fs>
   int fieldbuffer_get ( FieldBufferPyObject * exporter, Py_buffer * view, int flags )
   {
      return FieldBufferGetDispatch<T,fs>::get( exporter, view, flags, false );
   }


   template<typename T, uint_t fs>
   int fieldbuffer_release ( PyObject * /*exporter*/, Py_buffer * view )
   {
      delete [] view->strides;
      delete [] view->shape;
      //std::cout << "Releasing Field Buffer " << std::endl;
      return 0;
   }


   template<typename T, uint_t fs>
   static void fieldbuffer_dealloc( FieldBufferPyObject * exporter )
   {
      namespace bp = boost::python;
      bp::object fieldAllocObject ( bp::handle<>(bp::borrowed( exporter->fieldAlloc ) ) );
      FieldAllocator<T> * fieldAlloc = bp::extract< FieldAllocator<T> * > ( fieldAllocObject );

      fieldAlloc->decrementReferenceCount( (T*) exporter->fieldData );

      Py_DECREF( exporter->fieldAlloc );


      //std::cout << "Dealloc Fieldbuffer " << (void*) exporter << std::endl;
      Py_TYPE(exporter)->tp_free ((PyObject*) exporter );
   }


   template<typename T, uint_t fs>
   Py_ssize_t fieldbuffer_getbuffer(FieldBufferPyObject *, Py_ssize_t , const void **)
   {
      WALBERLA_CHECK(false, "fieldbuffer_getbuffer is part of the old buffer interface and should never be used");
      return Py_ssize_t(0);// prevent compiler warning
   }


   template<typename T, uint_t fs>
   Py_ssize_t fieldbuffer_getsegcount(FieldBufferPyObject *, Py_ssize_t *)
   {
      WALBERLA_CHECK(false, "fieldbuffer_getsegcount is part of the old buffer interface and should never be used");
      return Py_ssize_t(0);// prevent compiler warning
   }


   template<typename T, uint_t fs>
   Py_ssize_t fieldbuffer_getbuffer_withGl(FieldBufferPyObject *, Py_ssize_t , const void **)
   {
      WALBERLA_CHECK(false, "fieldbuffer_getbuffer is part of the old buffer interface and should never be used");
      return Py_ssize_t(0);// prevent compiler warning
   }


   template<typename T, uint_t fs>
   Py_ssize_t fieldbuffer_getsegcount_withGl(FieldBufferPyObject *, Py_ssize_t *)
   {
      WALBERLA_CHECK(false, "fieldbuffer_getsegcount is part of the old buffer interface and should never be used");
      return Py_ssize_t(0);// prevent compiler warning
   }



   #ifdef WALBERLA_CXX_COMPILER_IS_GNU
   #pragma GCC diagnostic push
   #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
   #endif


   template<typename T, uint_t fs>
   struct FieldBufferType
   {
      static PyBufferProcs bufferProcs;
      static PyTypeObject  value;
   };

   template<typename T, uint_t fs>
   PyBufferProcs FieldBufferType<T,fs>::bufferProcs = {
#if PY_MAJOR_VERSION < 3
       (readbufferproc)   fieldbuffer_getbuffer  <T,fs>, /* bf_getreadbuffer */
       (writebufferproc)  fieldbuffer_getbuffer  <T,fs>, /* bf_getwritebuffer */
       (segcountproc)     fieldbuffer_getsegcount<T,fs>, /* bf_getsegcount */
       (charbufferproc)   fieldbuffer_getbuffer  <T,fs>,  /* bf_getcharbuffer */
#endif
       (getbufferproc)    fieldbuffer_get    <T,fs>, /* bf_getbuffer */
       (releasebufferproc)fieldbuffer_release<T,fs>, /* bf_releasebuffer */
   };

   template<typename T, uint_t fs>
   PyTypeObject FieldBufferType<T,fs>::value = {
       PyVarObject_HEAD_INIT(NULL, 0)
       "walberla_cpp.FieldBuffer",              /* tp_name */
       sizeof(FieldBufferPyObject),             /* tp_basicsize */
       0,                                       /* tp_itemsize */
       (destructor)fieldbuffer_dealloc<T,fs>,   /* tp_dealloc */
       0,                                       /* tp_print */
       0,                                       /* tp_getattr */
       0,                                       /* tp_setattr */
       0,                                       /* tp_reserved */
       0,                                       /* tp_repr */
       0,                                       /* tp_as_number */
       0,                                       /* tp_as_sequence */
       0,                                       /* tp_as_mapping */
       0,                                       /* tp_hash  */
       0,                                       /* tp_call */
       0,                                       /* tp_str */
       0,                                       /* tp_getattro */
       0,                                       /* tp_setattro */
       &FieldBufferType<T,fs>::bufferProcs,     /* tp_as_buffer */
#if PY_MAJOR_VERSION >= 3
       Py_TPFLAGS_DEFAULT,                      /* tp_flags */
#else
       Py_TPFLAGS_DEFAULT |
          Py_TPFLAGS_HAVE_NEWBUFFER,            /* tp_flags */
#endif
       "FieldBuffer Objects",                   /* tp_doc */
   };



   template<typename T, uint_t fs>
   struct FieldBufferTypeGl
   {
      static PyBufferProcs bufferProcs;
      static PyTypeObject  value;
   };

   template<typename T, uint_t fs>
   PyBufferProcs FieldBufferTypeGl<T,fs>::bufferProcs = {
#if PY_MAJOR_VERSION < 3
       (readbufferproc)   fieldbuffer_getbuffer_withGl  <T,fs>, /* bf_getreadbuffer */
       (writebufferproc)  fieldbuffer_getbuffer_withGl  <T,fs>, /* bf_getwritebuffer */
       (segcountproc)     fieldbuffer_getsegcount_withGl<T,fs>, /* bf_getsegcount */
       (charbufferproc)   fieldbuffer_getbuffer_withGl  <T,fs>,  /* bf_getcharbuffer */
#endif
       (getbufferproc)    fieldbuffer_get_withGl<T,fs>, /* bf_getbuffer */
       (releasebufferproc)fieldbuffer_release   <T,fs>, /* bf_releasebuffer */
   };

   template<typename T, uint_t fs>
   PyTypeObject FieldBufferTypeGl<T,fs>::value = {
       PyVarObject_HEAD_INIT(NULL, 0)
       "walberla_cpp.FieldBufferGl",            /* tp_name */
       sizeof(FieldBufferPyObject),             /* tp_basicsize */
       0,                                       /* tp_itemsize */
       (destructor)fieldbuffer_dealloc<T,fs>,   /* tp_dealloc */
       0,                                       /* tp_print */
       0,                                       /* tp_getattr */
       0,                                       /* tp_setattr */
       0,                                       /* tp_reserved */
       0,                                       /* tp_repr */
       0,                                       /* tp_as_number */
       0,                                       /* tp_as_sequence */
       0,                                       /* tp_as_mapping */
       0,                                       /* tp_hash  */
       0,                                       /* tp_call */
       0,                                       /* tp_str */
       0,                                       /* tp_getattro */
       0,                                       /* tp_setattro */
       &FieldBufferTypeGl<T,fs>::bufferProcs,   /* tp_as_buffer */
#if PY_MAJOR_VERSION >= 3
       Py_TPFLAGS_DEFAULT,                      /* tp_flags */
#else
       Py_TPFLAGS_DEFAULT |
          Py_TPFLAGS_HAVE_NEWBUFFER,            /* tp_flags */
#endif
       "FieldBufferGl Objects",                 /* tp_doc */
   };



   #ifdef WALBERLA_CXX_COMPILER_IS_GNU
   #pragma GCC diagnostic pop
   #endif

   // this will become the buffer() function of the field which creates the
   // FieldBuffer object, so this function is between field (exported in boost::python) and
   // FieldBuffer which is exported in native Python C API
   template<typename T, uint_t fs>
   boost::python::object field_getBufferInterface( boost::python::object field, bool withGhostLayers )
   {
      namespace bp = boost::python;

      FieldBufferPyObject *obj;
      if ( withGhostLayers )
         obj = (FieldBufferPyObject*) PyObject_CallObject((PyObject *) & FieldBufferTypeGl<T,fs>::value, NULL );
      else
         obj = (FieldBufferPyObject*) PyObject_CallObject((PyObject *) & FieldBufferType<T,fs>::value, NULL );


      Field<T,fs> * fieldPtr = bp::extract< Field<T,fs> * > ( field );
      bp::object fieldPtrObject( fieldPtr->getAllocator() );

      obj->field      = field.ptr();
      obj->fieldAlloc = fieldPtrObject.ptr();
      obj->fieldData  = (void*) ( fieldPtr->data() );
      Py_INCREF( obj->field      );
      Py_INCREF( obj->fieldAlloc );

      return bp::object ( bp::handle<>( (PyObject*) obj ) );
   }


   //===================================================================================================================
   //
   //  Aligned Allocation
   //
   //===================================================================================================================

   template<typename T>
   shared_ptr<field::FieldAllocator<T> > getAllocator(uint_t alignment)
   {
      if( alignment == 0 )
         return shared_ptr<field::FieldAllocator<T> >(); // leave to default - auto-detection of alignment
      else if ( alignment == 16 )
         return make_shared< field::AllocateAligned<T, 16> >();
      else if ( alignment == 32 )
         return make_shared< field::AllocateAligned<T, 32> >();
      else if ( alignment == 64 )
         return make_shared< field::AllocateAligned<T, 64> >();
      else if ( alignment == 128 )
         return make_shared< field::AllocateAligned<T, 128> >();
      else {
         PyErr_SetString( PyExc_ValueError, "Alignment parameter has to be one of 0, 16, 32, 64, 128." );
         throw boost::python::error_already_set();
         return shared_ptr<field::FieldAllocator<T> >();
      }
   }

   template< typename GhostLayerField_T >
   class GhostLayerFieldDataHandling : public field::BlockDataHandling< GhostLayerField_T >
   {
   public:
      typedef typename GhostLayerField_T::value_type Value_T;

      GhostLayerFieldDataHandling( const weak_ptr<StructuredBlockStorage> &blocks, const uint_t nrOfGhostLayers,
                                   const Value_T &initValue, const Layout layout, uint_t alignment = 0 ) :
              blocks_( blocks ), nrOfGhostLayers_( nrOfGhostLayers ), initValue_( initValue ), layout_( layout ),
              alignment_( alignment ) {}

      GhostLayerField_T * allocate( IBlock * const block )
      {
         auto blocks = blocks_.lock();
         WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'AlwaysInitializeBlockDataHandling' for a block "
                                             "storage object that doesn't exist anymore" );
         GhostLayerField_T * field = new GhostLayerField_T ( blocks->getNumberOfXCells( *block ),
                                                             blocks->getNumberOfYCells( *block ),
                                                             blocks->getNumberOfZCells( *block ),
                                                             nrOfGhostLayers_, initValue_, layout_,
                                                             getAllocator<Value_T>(alignment_) );
         return field;
      }

      GhostLayerField_T * reallocate( IBlock * const block )
      {
         return allocate(block);
      }

   private:
      weak_ptr< StructuredBlockStorage > blocks_;

      uint_t  nrOfGhostLayers_;
      Value_T initValue_;
      Layout  layout_;
      uint_t alignment_;
   };


   //===================================================================================================================
   //
   //  Field functions redefined for easier export
   //
   //===================================================================================================================


   static inline Cell tupleToCell( boost::python::tuple & tuple  )
   {
      using boost::python::extract;
      return Cell (  extract<cell_idx_t>( tuple[0] ),
                     extract<cell_idx_t>( tuple[1] ),
                     extract<cell_idx_t>( tuple[2] ) );
   }

   template<typename Field_T>
   void field_setCellXYZ( Field_T & field, boost::python::tuple args, const typename Field_T::value_type & value )
   {
      using namespace boost::python;

      if ( len(args) < 3 || len(args) > 4 )
      {
         PyErr_SetString( PyExc_RuntimeError, "3 or 4 indices required");
         throw error_already_set();
      }

      cell_idx_t f = 0;
      if ( len(args) == 4 )
         f = extract<cell_idx_t> ( args[3] );

      Cell cell = tupleToCell(args);
      if ( ! field.coordinatesValid( cell[0], cell[1], cell[2], f ) )
      {
         PyErr_SetString( PyExc_IndexError, "Field indices out of bounds");
         throw error_already_set();
      }
      field(cell, f) = value;
   }

   template<typename Field_T>
   typename Field_T::value_type field_getCellXYZ( Field_T & field, boost::python::tuple args )
   {
      using namespace boost::python;
      if ( len(args) < 3 || len(args) > 4 )
      {
         PyErr_SetString( PyExc_RuntimeError, "3 or 4 indices required");
         throw error_already_set();
      }

      cell_idx_t f = 0;
      if ( len(args) == 4 )
         f = extract<cell_idx_t> ( args[3] );

      Cell cell = tupleToCell(args);
      if ( ! field.coordinatesValid( cell[0], cell[1], cell[2], f ) )
      {
         PyErr_SetString( PyExc_IndexError, "Field indices out of bounds");
         throw error_already_set();
      }

      return field( cell, f );
   }


   template<typename Field_T>
   boost::python::object field_size( const Field_T & field ) {
      return boost::python::make_tuple( field.xSize(), field.ySize(), field.zSize(), field.fSize() );
   }

   template<typename GlField_T>
   boost::python::object field_sizeWithGhostLayer( const GlField_T & field ) {
      return boost::python::make_tuple( field.xSizeWithGhostLayer(), field.ySizeWithGhostLayer(), field.zSizeWithGhostLayer(), field.fSize() );
   }


   template<typename Field_T>
   boost::python::object field_allocSize( const Field_T & field ) {
      return boost::python::make_tuple( field.xAllocSize(), field.yAllocSize(),
                                        field.zAllocSize(), field.fAllocSize() );
   }


   template<typename Field_T>
   boost::python::object field_strides( const Field_T & field ) {
      return boost::python::make_tuple( field.xStride(), field.yStride(),
                                        field.zStride(), field.fStride() );
   }

   template<typename Field_T>
   boost::python::object field_offsets( const Field_T & field ) {
      return boost::python::make_tuple( field.xOff(), field.yOff(), field.zOff() );
   }


   template<typename Field_T>
   boost::python::object field_layout( const Field_T & f ) {
      if ( f.layout() == field::fzyx ) return boost::python::object( "fzyx" );
      if ( f.layout() == field::zyxf ) return boost::python::object( "zyxf" );

      return boost::python::object();
   }


   template<typename Field_T>
   void field_swapDataPointers( Field_T & f1, Field_T & f2 )
   {
      if ( ! f1.hasSameAllocSize(f2 ) ||
           ! f1.hasSameSize( f2)      ||
             f1.layout() != f2.layout() )
      {
         PyErr_SetString( PyExc_ValueError, "The data of fields with different sizes or layout cannot be swapped");
         throw boost::python::error_already_set();
      }
      f1.swapDataPointers( f2 );
   }


   template<typename T>
   T  FF_getFlag( const FlagField<T> & ff,  const std::string & flag ) {
      if ( ! ff.flagExists(flag) )
      {
         PyErr_SetString( PyExc_ValueError, "No such flag");
         throw boost::python::error_already_set();
      }
      return ff.getFlag( flag );
   }

   template<typename T>
   boost::python::object FF_registeredFlags( const FlagField<T> & ff )
   {
      std::vector<FlagUID> flags;
      ff.getAllRegisteredFlags( flags );
      boost::python::list result;

      for( auto i = flags.begin(); i != flags.end(); ++i )
         result.append( i->toString() );
      boost::python::object objectResult = result;
      return objectResult;
   }


   template<typename T>
   boost::python::object FF_flagMap( const FlagField<T> & ff )
   {
      std::vector<FlagUID> flags;
      ff.getAllRegisteredFlags( flags );
      boost::python::dict result;

      for( auto i = flags.begin(); i != flags.end(); ++i )
         result[ i->toString() ] = ff.getFlag( *i );
      boost::python::object objectResult = result;
      return objectResult;
   }

   template<typename T>
   boost::python::object FF_registerFlag( FlagField<T> & ff, const std::string & flag, boost::python::object bitNr )
   {
      using namespace boost::python;

      try {

         if ( bitNr == object() )
            return object( ff.registerFlag( FlagUID(flag) ) );
         else
         {
            if ( extract<uint_t>(bitNr ).check() ) {
               uint_t bit = extract<uint_t>(bitNr);
               return object( ff.registerFlag( flag, bit ) );
            }
            else {
               PyErr_SetString( PyExc_ValueError, "Parameter bitNr has to be a positive integer");
               throw boost::python::error_already_set();
            }
         }
      }
      catch ( std::runtime_error & e ) {
         PyErr_SetString( PyExc_ValueError, e.what() );
         throw boost::python::error_already_set();
      }

   }

   template<typename T>
   std::string FF_getFlagName( const FlagField<T> & ff, T flag )
   {
      try {
         return ff.getFlagUID( flag ).getIdentifier();
      }
      catch ( std::runtime_error & e ) {
         PyErr_SetString( PyExc_ValueError, e.what() );
         throw boost::python::error_already_set();
      }
   }


   template<typename Field_T>
   boost::python::object copyAdaptorToField( const Field_T & f )
   {
      typedef GhostLayerField<typename Field_T::value_type, Field_T::F_SIZE> ResField;
      auto res = make_shared< ResField > ( f.xSize(), f.ySize(), f.zSize(), f.nrOfGhostLayers() );

      auto srcIt = f.beginWithGhostLayerXYZ();
      auto dstIt = res->beginWithGhostLayerXYZ();
      while ( srcIt != f.end() )
      {
         for( cell_idx_t fCoord = 0; fCoord < cell_idx_c(Field_T::F_SIZE); ++fCoord )
            dstIt.getF( fCoord ) = srcIt.getF( fCoord );

         ++srcIt;
         ++dstIt;
      }
      return boost::python::object( res );
   }



   //===================================================================================================================
   //
   //  Field export
   //
   //===================================================================================================================

   template<typename T>
   void exportFlagFieldIfUnsigned( typename std::enable_if<std::is_unsigned<T>::value >::type* = 0 )
   {
      using namespace boost::python;

      class_< FlagField<T> ,
              shared_ptr<FlagField<T> >,
              bases<GhostLayerField<T,1> >,
              boost::noncopyable > ( "FlagField", no_init )
          .def         ( "registerFlag", &FF_registerFlag<T>   ,  ( arg("flagName"),  arg("bitNr") = object() ) )
          .def         ( "flag",         &FF_getFlag<T>         )
          .def         ( "flagName",     &FF_getFlagName<T>     )
          .add_property( "flags",        &FF_registeredFlags<T> )
          .add_property( "flagMap",      &FF_flagMap<T>         )
          ;

   }
   template<typename T>
   void exportFlagFieldIfUnsigned( typename std::enable_if< ! std::is_unsigned<T>::value >::type* = 0 )  {}


   struct FieldExporter
   {
      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         typedef typename FieldType::value_type T;
         const uint_t F_SIZE = FieldType::F_SIZE;
         typedef GhostLayerField<T,F_SIZE> GlField_T;
         typedef Field<T,F_SIZE> Field_T;

         using namespace boost::python;

         class_<Field_T, shared_ptr<Field_T>, boost::noncopyable>( "Field", no_init )
            .add_property("layout",    &field_layout          < Field_T > )
            .add_property("size",      &field_size            < Field_T > )
            .add_property("allocSize", &field_allocSize       < Field_T > )
            .add_property("strides",   &field_strides         < Field_T > )
            .add_property("offsets",   &field_offsets         < Field_T > )
            .def("clone",              &Field_T::clone             , return_value_policy<manage_new_object>())
            .def("cloneUninitialized", &Field_T::cloneUninitialized, return_value_policy<manage_new_object>())
            .def("swapDataPointers",   &field_swapDataPointers< Field_T > )
            .def("__getitem__",        &field_getCellXYZ      < Field_T > )
            .def("__setitem__",        &field_setCellXYZ      < Field_T > )
            .def("buffer",             &field_getBufferInterface<T,F_SIZE>, ( arg("withGhostLayers") = false ) );
            ;

         class_< GlField_T , shared_ptr<GlField_T>, bases<Field_T>,  boost::noncopyable > ( "GhostLayerField", no_init )
            .add_property("sizeWithGhostLayer", &field_sizeWithGhostLayer<GlField_T> )
            .add_property("nrOfGhostLayers",    &GlField_T::nrOfGhostLayers     )
         ;

         if ( F_SIZE == 1 )
            exportFlagFieldIfUnsigned<T>();


         // Field Buffer
         FieldBufferType<T,F_SIZE>::value.tp_new = PyType_GenericNew;
         if ( PyType_Ready(& FieldBufferType<T,F_SIZE>::value ) < 0 )
             return;

         Py_INCREF( (& FieldBufferType<T,F_SIZE>::value ) );
         PyModule_AddObject( boost::python::scope().ptr(), "FieldBuffer", (PyObject *)&FieldBufferType<T,F_SIZE>::value );

         // Field Buffer with ghost layer
         FieldBufferTypeGl<T,F_SIZE>::value.tp_new = PyType_GenericNew;
         if ( PyType_Ready(& FieldBufferTypeGl<T,F_SIZE>::value ) < 0 )
             return;

         Py_INCREF( (& FieldBufferTypeGl<T,F_SIZE>::value ) );
         PyModule_AddObject( boost::python::scope().ptr(), "FieldBufferGl", (PyObject *)&FieldBufferTypeGl<T,F_SIZE>::value );


         // Field Buffer type

         if (  python_coupling::isTypeRegisteredInBoostPython< FieldAllocator<T> >() == false )
         {
            class_< FieldAllocator<T>, shared_ptr<FieldAllocator<T> >, boost::noncopyable> ( "FieldAllocator", no_init )
                  .def( "incrementReferenceCount", &FieldAllocator<T>::incrementReferenceCount  )
                  .def( "decrementReferenceCount", &FieldAllocator<T>::decrementReferenceCount  );
         }

         using field::communication::PackInfo;
         class_< PackInfo<GlField_T>,
                 shared_ptr< PackInfo<GlField_T> >,
                 bases<walberla::communication::UniformPackInfo>,
                 boost::noncopyable >( "FieldPackInfo", no_init );


         using field::communication::UniformMPIDatatypeInfo;
         class_< UniformMPIDatatypeInfo<GlField_T>,
                 shared_ptr< UniformMPIDatatypeInfo<GlField_T> >,
                 bases<walberla::communication::UniformMPIDatatypeInfo>,
                 boost::noncopyable >( "FieldMPIDataTypeInfo", no_init );

      }
   };


   struct GhostLayerFieldAdaptorExporter
   {
      GhostLayerFieldAdaptorExporter(  const std::string & name )
         : name_ ( name )
      {}

      template< typename Adaptor>
      void operator()( python_coupling::NonCopyableWrap<Adaptor> )
      {
         using namespace boost::python;

         class_< Adaptor, shared_ptr<Adaptor>, boost::noncopyable> ( name_.c_str(), no_init )
               .add_property("size",                &field_size              <Adaptor > )
               .add_property("sizeWithGhostLayer",  &field_sizeWithGhostLayer<Adaptor> )
               .add_property("nrOfGhostLayers",     &Adaptor::nrOfGhostLayers          )
               .def("__getitem__",                  &field_getCellXYZ        <Adaptor> )
               .def("copyToField",                  &copyAdaptorToField      <Adaptor> );
      }

      std::string name_;
   };

   //===================================================================================================================
   //
   //  createField
   //
   //===================================================================================================================


   class CreateFieldExporter
   {
   public:
      CreateFieldExporter( uint_t xs, uint_t ys, uint_t zs, uint_t fs, uint_t gl,
                           Layout layout, const boost::python::object & type, uint_t alignment,
                           const shared_ptr<boost::python::object> & resultPointer  )
         : xs_( xs ), ys_(ys), zs_(zs), fs_(fs), gl_(gl),
           layout_( layout),  type_( type ), alignment_(alignment), resultPointer_( resultPointer )
      {}

      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         using namespace boost::python;
         typedef typename FieldType::value_type T;
         const uint_t F_SIZE = FieldType::F_SIZE;

         if( F_SIZE != fs_ )
            return;

         if( python_coupling::isCppEqualToPythonType<T>( (PyTypeObject *)type_.ptr() )  )
         {
            T initVal = T(); //extract<T> ( initValue_ );
            *resultPointer_ = object( make_shared< GhostLayerField<T,F_SIZE> >( xs_,ys_,zs_, gl_, initVal, layout_,
                                                                                getAllocator<T>(alignment_)));
         }
      }

   private:
      uint_t xs_;
      uint_t ys_;
      uint_t zs_;
      uint_t fs_;
      uint_t gl_;
      Layout layout_;
      boost::python::object type_;
      uint_t alignment_;
      shared_ptr<boost::python::object> resultPointer_;
   };

   template<typename FieldTypes>
   boost::python::object createPythonField( boost::python::list size,
                                            boost::python::object type,
                                            uint_t ghostLayers,
                                            Layout layout,
                                            uint_t alignment)
   {
      using namespace boost::python;
      uint_t xSize = extract<uint_t> ( size[0] );
      uint_t ySize = extract<uint_t> ( size[1] );
      uint_t zSize = extract<uint_t> ( size[2] );
      uint_t sizeLen = uint_c( len( size ) );
      uint_t fSize = 1;
      if ( sizeLen == 4 )
         fSize = extract<uint_t> ( size[3] );

      if ( ! PyType_Check( type.ptr() ) ) {
         PyErr_SetString( PyExc_RuntimeError, "Invalid 'type' parameter");
         throw error_already_set();
      }

      auto result = make_shared<boost::python::object>();
      CreateFieldExporter exporter( xSize,ySize, zSize, fSize, ghostLayers, layout, type, alignment, result );
      python_coupling::for_each_noncopyable_type< FieldTypes >  ( exporter );

      if ( *result == object()  )
      {
         PyErr_SetString( PyExc_ValueError, "Cannot create field of this (type,f-size) combination");
         throw error_already_set();
      }
      else {
         return *result;
      }
   }

   //===================================================================================================================
   //
   //  createFlagField
   //
   //===================================================================================================================


   inline boost::python::object createPythonFlagField( boost::python::list size, uint_t nrOfBits, uint_t ghostLayers )
   {
      using namespace boost::python;

      uint_t sizeLen = uint_c( len( size ) );
      if ( sizeLen != 3 ) {
         PyErr_SetString( PyExc_ValueError, "Size parameter has to be a list of length 3");
         throw error_already_set();
      }
      uint_t xSize = extract<uint_t> ( size[0] );
      uint_t ySize = extract<uint_t> ( size[1] );
      uint_t zSize = extract<uint_t> ( size[2] );

      if( nrOfBits == 8 )        return object( make_shared< FlagField< uint8_t > >( xSize, ySize, zSize, ghostLayers )  );
      else if ( nrOfBits == 16 ) return object( make_shared< FlagField< uint16_t> >( xSize, ySize, zSize, ghostLayers )  );
      else if ( nrOfBits == 32 ) return object( make_shared< FlagField< uint32_t> >( xSize, ySize, zSize, ghostLayers )  );
      else if ( nrOfBits == 64 ) return object( make_shared< FlagField< uint64_t> >( xSize, ySize, zSize, ghostLayers )  );
      else
      {
         PyErr_SetString( PyExc_ValueError, "Allowed values for number of bits are: 8,16,32,64");
         throw error_already_set();
      }
   }


   //===================================================================================================================
   //
   //  addToStorage
   //
   //===================================================================================================================

   class AddToStorageExporter
   {
   public:
      AddToStorageExporter(const shared_ptr<StructuredBlockStorage> & blocks,
                           const std::string & name, uint_t fs, uint_t gl, Layout layout,
                           const boost::python::object & type,
                           const boost::python::object & initObj,
                           uint_t alignment )
         : blocks_( blocks ), name_( name ), fs_( fs ),
           gl_(gl),layout_( layout),  type_( type ), initObj_( initObj), alignment_(alignment), found_(false)
      {}

      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         using namespace boost::python;
         typedef typename FieldType::value_type T;
         const uint_t F_SIZE = FieldType::F_SIZE;

         if( F_SIZE != fs_ )
            return;

         if( !found_ && python_coupling::isCppEqualToPythonType<T>( (PyTypeObject *)type_.ptr() )  )
         {
            typedef internal::GhostLayerFieldDataHandling< GhostLayerField<T,F_SIZE > > DataHandling;
            if ( initObj_ == object() ) {
               auto dataHandling = walberla::make_shared< DataHandling >( blocks_, gl_, T(), layout_, alignment_ );
               blocks_->addBlockData( dataHandling, name_ );
            }
            else {
               auto dataHandling = walberla::make_shared< DataHandling >( blocks_, gl_, extract<T>(initObj_), layout_, alignment_ );
               blocks_->addBlockData( dataHandling, name_ );
            }
            found_ = true;
         }
      }

      bool successful() const { return found_; }
   private:
      shared_ptr< StructuredBlockStorage > blocks_;
      std::string name_;
      uint_t fs_;
      uint_t gl_;
      Layout layout_;
      boost::python::object type_;
      boost::python::object initObj_;
      uint_t alignment_;
      bool found_;
   };

   template<typename FieldTypes>
   void addToStorage( const shared_ptr<StructuredBlockStorage> & blocks, const std::string & name,
                      boost::python::object type, uint_t fs, uint_t gl, Layout layout, boost::python::object initValue,
                      uint_t alignment)
   {
      using namespace boost::python;

      if ( ! PyType_Check( type.ptr() ) ) {
         PyErr_SetString( PyExc_RuntimeError, "Invalid 'type' parameter");
         throw error_already_set();
      }

      auto result = make_shared<boost::python::object>();
      AddToStorageExporter exporter( blocks, name, fs, gl, layout, type, initValue, alignment );
      python_coupling::for_each_noncopyable_type< FieldTypes >  ( std::ref(exporter) );

      if ( ! exporter.successful() ) {
         PyErr_SetString( PyExc_ValueError, "Adding Field failed.");
         throw error_already_set();
      }
   }


   inline void addFlagFieldToStorage( const shared_ptr<StructuredBlockStorage> & blocks, const std::string & name,
                              uint_t nrOfBits, uint_t gl  )
   {
      if( nrOfBits == 8 )        field::addFlagFieldToStorage< FlagField<uint8_t>  > ( blocks, name, gl );
      else if ( nrOfBits == 16 ) field::addFlagFieldToStorage< FlagField<uint16_t> > ( blocks, name, gl );
      else if ( nrOfBits == 32 ) field::addFlagFieldToStorage< FlagField<uint32_t> > ( blocks, name, gl );
      else if ( nrOfBits == 64 ) field::addFlagFieldToStorage< FlagField<uint64_t> > ( blocks, name, gl );
      else
      {
         PyErr_SetString( PyExc_ValueError, "Allowed values for number of bits are: 8,16,32,64");
         throw boost::python::error_already_set();
      }
   }

   //===================================================================================================================
   //
   //  createVTKWriter
   //
   //===================================================================================================================

   class CreateVTKWriterExporter
   {
   public:
       CreateVTKWriterExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                                ConstBlockDataID fieldId, const std::string & vtkName)
         : blocks_( blocks ), fieldId_(fieldId), vtkName_( vtkName )
      {}

      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         using namespace boost::python;

         IBlock * firstBlock =  & ( * blocks_->begin() );
         if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
            writer_ = shared_ptr<field::VTKWriter<FieldType> >( new field::VTKWriter<FieldType>(fieldId_, vtkName_));
      }

      shared_ptr< vtk::BlockCellDataWriterInterface > getCreatedWriter() {
         return writer_;
      }

   private:
      shared_ptr< vtk::BlockCellDataWriterInterface > writer_;
      shared_ptr< StructuredBlockStorage > blocks_;
      ConstBlockDataID fieldId_;
      std::string vtkName_;
   };


   template<typename FieldTypes>
   inline shared_ptr<vtk::BlockCellDataWriterInterface> createVTKWriter(const shared_ptr<StructuredBlockStorage> & blocks,
                                                                        const std::string & name,
                                                                        const std::string & nameInVtkOutput = "")
   {
      std::string vtkName = nameInVtkOutput;
      if( vtkName.size() == 0)
         vtkName = name;

      if ( blocks->begin() == blocks->end() )
         return shared_ptr<vtk::BlockCellDataWriterInterface>();
      auto fieldID = python_coupling::blockDataIDFromString( *blocks, name );

      CreateVTKWriterExporter exporter(blocks, fieldID, vtkName);
      python_coupling::for_each_noncopyable_type< FieldTypes >  ( std::ref(exporter) );
      if ( ! exporter.getCreatedWriter() ) {
         PyErr_SetString( PyExc_ValueError, "Failed to create writer");
         throw boost::python::error_already_set();
      }
      else {
         return exporter.getCreatedWriter();
      }
   }

   //===================================================================================================================
   //
   //  createFlagFieldVTKWriter
   //
   //===================================================================================================================

   class CreateFlagFieldVTKWriterExporter
   {
   public:
      CreateFlagFieldVTKWriterExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                                        ConstBlockDataID fieldId, const std::string & vtkName,
                                        boost::python::dict flagMapping)
         : blocks_( blocks ), fieldId_(fieldId), vtkName_( vtkName ), flagMapping_( flagMapping )
      {}

      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         using namespace boost::python;

         IBlock * firstBlock =  & ( * blocks_->begin() );
         if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
         {
            typedef typename FieldType::flag_t flag_t;
            typedef field::FlagFieldMapping<FieldType, flag_t> FFMapping;
            auto uncastedWriter = shared_ptr<FFMapping >( new FFMapping(fieldId_, vtkName_));
            writer_ = uncastedWriter;
            auto keys = flagMapping_.keys();

            for( int i=0; i < len(keys); ++i ) {
               uncastedWriter->addMapping(FlagUID(extract<std::string>(keys[i])),
                                          extract<flag_t>(flagMapping_[keys[i]]) );
            }
         }

      }

      shared_ptr< vtk::BlockCellDataWriterInterface > getCreatedWriter() {
         return writer_;
      }

   private:
      shared_ptr< vtk::BlockCellDataWriterInterface > writer_;
      shared_ptr< StructuredBlockStorage > blocks_;
      ConstBlockDataID fieldId_;
      std::string vtkName_;
      boost::python::dict flagMapping_;
   };


   template<typename FieldTypes>
   inline shared_ptr<vtk::BlockCellDataWriterInterface> createFlagFieldVTKWriter(const shared_ptr<StructuredBlockStorage> & blocks,
                                                                                 const std::string & name,
                                                                                 boost::python::dict flagMapping,
                                                                                 const std::string & nameInVtkOutput = "" )
   {
      std::string vtkName = nameInVtkOutput;
      if( vtkName.size() == 0)
         vtkName = name;

      if ( blocks->begin() == blocks->end() )
         return shared_ptr<vtk::BlockCellDataWriterInterface>();
      auto fieldID = python_coupling::blockDataIDFromString( *blocks, name );

      CreateFlagFieldVTKWriterExporter exporter(blocks, fieldID, vtkName, flagMapping);
      python_coupling::for_each_noncopyable_type< FieldTypes >  ( std::ref(exporter) );
      if ( ! exporter.getCreatedWriter() ) {
         PyErr_SetString( PyExc_ValueError, "Failed to create writer");
         throw boost::python::error_already_set();
      }
      else {
         return exporter.getCreatedWriter();
      }
   }


   //===================================================================================================================
   //
   //  createBinarizationFieldWriter
   //
   //===================================================================================================================

   class CreateBinarizationVTKWriterExporter
   {
   public:
      CreateBinarizationVTKWriterExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                                           ConstBlockDataID fieldId, const std::string & vtkName, uint_t mask )
         : blocks_( blocks ), fieldId_(fieldId), vtkName_( vtkName ), mask_(mask)
      {}

      template< typename FieldType>
      void operator() ( python_coupling::NonCopyableWrap<FieldType> )
      {
         using namespace boost::python;

         IBlock * firstBlock =  & ( * blocks_->begin() );
         if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
         {
            typedef field::BinarizationFieldWriter<FieldType> Writer;
            writer_ = shared_ptr<Writer>( new Writer(fieldId_, vtkName_,
                                                     static_cast<typename FieldType::value_type>(mask_) ));
         }
      }

      shared_ptr< vtk::BlockCellDataWriterInterface > getCreatedWriter() {
         return writer_;
      }

   private:
      shared_ptr< vtk::BlockCellDataWriterInterface > writer_;
      shared_ptr< StructuredBlockStorage > blocks_;
      ConstBlockDataID fieldId_;
      std::string vtkName_;
      uint_t mask_;
   };


   template<typename FieldTypes>
   inline shared_ptr<vtk::BlockCellDataWriterInterface> createBinarizationVTKWriter(const shared_ptr<StructuredBlockStorage> & blocks,
                                                                                    const std::string & name,
                                                                                    uint_t mask,
                                                                                    const std::string & nameInVtkOutput = "" )
   {
      std::string vtkName = nameInVtkOutput;
      if( vtkName.size() == 0)
         vtkName = name;

      if ( blocks->begin() == blocks->end() )
         return shared_ptr<vtk::BlockCellDataWriterInterface>();
      auto fieldID = python_coupling::blockDataIDFromString( *blocks, name );

      CreateBinarizationVTKWriterExporter exporter(blocks, fieldID, vtkName, mask);
      python_coupling::for_each_noncopyable_type< FieldTypes >  ( std::ref(exporter) );
      if ( ! exporter.getCreatedWriter() ) {
         PyErr_SetString( PyExc_ValueError, "Failed to create writer");
         throw boost::python::error_already_set();
      }
      else {
         return exporter.getCreatedWriter();
      }
   }


} // namespace internal




template<typename FieldTypes >
void exportFields()
{
   using namespace boost::python;

   enum_<Layout>("Layout")
       .value("fzyx", fzyx)
       .value("zyxf", zyxf)
       .export_values();

   python_coupling::for_each_noncopyable_type< FieldTypes > ( internal::FieldExporter() );

   def( "createField", &internal::createPythonField<FieldTypes>, ( ( arg("size")                    ),
                                                                   ( arg("type")                    ),
                                                                   ( arg("ghostLayers") = uint_t(1) ),
                                                                   ( arg("layout")      = zyxf      ),
                                                                   ( arg("alignment")   = 0         )) );

   def( "createFlagField", &internal::createPythonFlagField, ( ( arg("size")                      ),
                                                               ( arg("nrOfBits")    = uint_t(32)  ),
                                                               ( arg("ghostLayers") = uint_t(1)   )  ) );

   def( "addToStorage",    &internal::addToStorage<FieldTypes>, ( ( arg("blocks")                  ),
                                                                  ( arg("name")                    ),
                                                                  ( arg("type")                    ),
                                                                  ( arg("fSize")       = 1         ),
                                                                  ( arg("ghostLayers") = uint_t(1) ),
                                                                  ( arg("layout")      = zyxf      ),
                                                                  ( arg("initValue")   = object()  ),
                                                                  ( arg("alignment")   = 0         ) ) );

   def( "addFlagFieldToStorage",&internal::addFlagFieldToStorage, ( ( arg("blocks")                  ),
                                                                    ( arg("name")                    ),
                                                                    ( arg("nrOfBits")=8              ),
                                                                    ( arg("ghostLayers") = uint_t(1) ) ) );

   def( "createVTKWriter", &internal::createVTKWriter<FieldTypes>, ( arg("blocks"), arg("name"), arg("vtkName")="" ));


   typedef boost::mpl::vector<
           FlagField<uint8_t>,
           FlagField<uint16_t>,
           FlagField<uint32_t>,
           FlagField<uint64_t> > FlagFields;

   def( "createFlagFieldVTKWriter", &internal::createFlagFieldVTKWriter<FlagFields>,
                                   ( arg("blocks"), arg("name"), arg("flagMapping"), arg("vtkName")="" ));


   typedef boost::mpl::vector<
           Field<uint8_t,1 >,
           Field<uint16_t, 1>,
           Field<uint32_t, 1>,
           Field<uint64_t, 1> > UintFields;

   def( "createBinarizationVTKWriter", &internal::createBinarizationVTKWriter<UintFields>,
        ( arg("blocks"), arg("name"), arg("mask"), arg("vtkName")="" ));
}


template<typename AdaptorTypes>
void exportGhostLayerFieldAdaptors()
{
   python_coupling::for_each_noncopyable_type< AdaptorTypes > ( internal::GhostLayerFieldAdaptorExporter("FieldAdaptor") );
}

template<typename AdaptorType>
void exportGhostLayerFieldAdaptor()
{
   typedef boost::mpl::vector< AdaptorType > AdaptorTypes;
   exportGhostLayerFieldAdaptors<AdaptorTypes>( );
}




} // namespace field
} // namespace walberla


