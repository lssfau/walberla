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
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#include "core/logging/Logging.h"
#include "cuda/GPUField.h"
#include "cuda/communication/GPUPackInfo.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/FieldCopy.h"
#include "cuda/GPUField.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "field/AddToStorage.h"
#include "field/python/FieldExport.h"
#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"

#include <type_traits>
#include <iostream>

namespace walberla {
namespace cuda {



namespace internal {

   //===================================================================================================================
   //
   //  Field export
   //
   //===================================================================================================================

   template<typename GpuField_T>
   uint64_t gpufield_ptr(const GpuField_T & gpuField)
   {
      return reinterpret_cast<uint64_t>(gpuField.pitchedPtr().ptr);
   }

   template<typename GpuField_T>
   std::string gpufield_dtypeStr(const GpuField_T & )
   {
      return std::string(field::internal::PythonFormatString<typename GpuField_T::value_type>::get());
   }

   struct GpuFieldExporter
   {
      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         using namespace boost::python;

         class_<GpuField_T, shared_ptr<GpuField_T>, boost::noncopyable>( "GpuField", no_init )
            .add_property("layout",              &field::internal::field_layout            < GpuField_T > )
            .add_property("size",                &field::internal::field_size              < GpuField_T > )
            .add_property("sizeWithGhostLayers", &field::internal::field_sizeWithGhostLayer< GpuField_T > )
            .add_property("allocSize",           &field::internal::field_allocSize         < GpuField_T > )
            .add_property("strides",             &field::internal::field_strides           < GpuField_T > )
            .add_property("offsets",             &field::internal::field_offsets           < GpuField_T > )
            .add_property("ptr",                 &gpufield_ptr                             < GpuField_T > )
            .add_property("dtypeStr",            &gpufield_dtypeStr                        < GpuField_T > )
            .add_property("isPitchedMem",        &GpuField_T::isPitchedMem )
            .def("swapDataPointers",             &field::internal::field_swapDataPointers  < GpuField_T > )
            .add_property("nrOfGhostLayers",     &GpuField_T::nrOfGhostLayers )
            .def("cloneUninitialized", &GpuField_T::cloneUninitialized, return_value_policy<manage_new_object>())
            ;


         using field::communication::PackInfo;
         using communication::GPUPackInfo;
         class_< GPUPackInfo<GpuField_T>,
                 shared_ptr< GPUPackInfo<GpuField_T> >,
                 bases<walberla::communication::UniformPackInfo>,
                 boost::noncopyable >( "GpuFieldPackInfo", no_init );

         using field::communication::UniformMPIDatatypeInfo;
         class_< UniformMPIDatatypeInfo<GpuField_T>,
                 shared_ptr< UniformMPIDatatypeInfo<GpuField_T> >,
                 bases<walberla::communication::UniformMPIDatatypeInfo>,
                 boost::noncopyable >( "GpuFieldMPIDataTypeInfo", no_init );

      }
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
                           Layout layout, const boost::python::object & type, bool usePitchedMem,
                           const shared_ptr<boost::python::object> & resultPointer )
         : xs_( xs ), ys_(ys), zs_(zs), fs_(fs), gl_(gl),
           layout_( layout),  type_( type ), usePitchedMem_( usePitchedMem ) , resultPointer_( resultPointer )
      {}

      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         using namespace boost::python;
         typedef typename GpuField_T::value_type T;
         if( python_coupling::isCppEqualToPythonType<T>( (PyTypeObject *)type_.ptr() )  )
         {
            *resultPointer_ = object( make_shared< GPUField<T> >( xs_,ys_,zs_, fs_,  gl_, layout_, usePitchedMem_ )  );
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
      bool usePitchedMem_;
      shared_ptr<boost::python::object> resultPointer_;
   };

   template<typename GpuFields>
   boost::python::object createPythonGpuField( boost::python::list size,
                                               boost::python::object type,
                                               uint_t ghostLayers,
                                               Layout layout,
                                               bool usePitchedMem)
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
      CreateFieldExporter exporter( xSize,ySize, zSize, fSize, ghostLayers, layout, type, usePitchedMem, result );
      python_coupling::for_each_noncopyable_type< GpuFields >( exporter );

      if ( *result == object()  )
      {
         PyErr_SetString( PyExc_ValueError, "Cannot create field of this type");
         throw error_already_set();
      }
      else {
         return *result;
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
      AddToStorageExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                           const std::string & name, uint_t fs, uint_t gl, Layout layout,
                           const boost::python::object & type,
                           bool usePitchedMem )
         : blocks_( blocks ), name_( name ), fs_( fs ),
           gl_(gl),layout_( layout),  type_( type ), usePitchedMem_(usePitchedMem), found_(false)
      {}

      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         typedef typename GpuField_T::value_type T;
         if( python_coupling::isCppEqualToPythonType<T>( (PyTypeObject *)type_.ptr() )  )
         {
            WALBERLA_ASSERT(!found_);
            addGPUFieldToStorage<GPUField<T> >(blocks_, name_, fs_, layout_, gl_, usePitchedMem_);
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
      bool usePitchedMem_;
      bool found_;
   };

   template<typename GpuFields>
   void addToStorage( const shared_ptr<StructuredBlockStorage> & blocks, const std::string & name,
                      boost::python::object type, uint_t fs, uint_t gl, Layout layout, bool usePitchedMem )
   {
      using namespace boost::python;

      if ( ! PyType_Check( type.ptr() ) ) {
         PyErr_SetString( PyExc_RuntimeError, "Invalid 'type' parameter");
         throw error_already_set();
      }

      auto result = make_shared<boost::python::object>();
      AddToStorageExporter exporter( blocks, name, fs, gl, layout, type, usePitchedMem );
      python_coupling::for_each_noncopyable_type<GpuFields>( std::ref(exporter) );

      if ( ! exporter.successful() ) {
         PyErr_SetString( PyExc_ValueError, "Adding Field failed.");
         throw error_already_set();
      }
   }


   //===================================================================================================================
   //
   //  createPackInfo Export
   //
   //===================================================================================================================

   template< typename GPUField_T >
   boost::python::object createGPUPackInfoToObject( BlockDataID bdId, uint_t numberOfGhostLayers )
   {
      using cuda::communication::GPUPackInfo;
      if ( numberOfGhostLayers > 0  )
         return boost::python::object( make_shared< GPUPackInfo<GPUField_T> >( bdId, numberOfGhostLayers ) );
      else
         return boost::python::object( make_shared< GPUPackInfo<GPUField_T> >( bdId ) );
   }

   FunctionExporterClass( createGPUPackInfoToObject, boost::python::object( BlockDataID, uint_t  ) );

   template< typename GpuFields>
   boost::python::object createPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                         const std::string & blockDataName, uint_t numberOfGhostLayers )
   {
      using cuda::communication::GPUPackInfo;

      auto bdId = python_coupling::blockDataIDFromString( *bs, blockDataName );
      if ( bs->begin() == bs->end() ) {
         // if no blocks are on this field an arbitrary PackInfo can be returned
         return createGPUPackInfoToObject< GPUField<real_t> > ( bdId, numberOfGhostLayers );
      }

      IBlock * firstBlock =  & ( * bs->begin() );
      python_coupling::Dispatcher<GpuFields, Exporter_createGPUPackInfoToObject > dispatcher( firstBlock );
      return dispatcher( bdId )( bdId, numberOfGhostLayers ) ;
   }


   //===================================================================================================================
   //
   //  createMPIDatatypeInfo
   //
   //===================================================================================================================


   template< typename GpuField_T >
   boost::python::object createMPIDatatypeInfoToObject( BlockDataID bdId, uint_t numberOfGhostLayers )
   {
      using field::communication::UniformMPIDatatypeInfo;
      if ( numberOfGhostLayers > 0 )
         return boost::python::object( make_shared< UniformMPIDatatypeInfo<GpuField_T> >( bdId, numberOfGhostLayers ) );
      else
         return boost::python::object( make_shared< UniformMPIDatatypeInfo<GpuField_T> >( bdId ) );
   }

   FunctionExporterClass( createMPIDatatypeInfoToObject, boost::python::object( BlockDataID, uint_t  ) );

   template< typename GpuFields>
   boost::python::object createMPIDatatypeInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                const std::string & blockDataName,
                                                uint_t numberOfGhostLayers)
   {
      auto bdId = python_coupling::blockDataIDFromString( *bs, blockDataName );
      if ( bs->begin() == bs->end() ) {
         // if no blocks are on this field an arbitrary MPIDatatypeInfo can be returned
         return createMPIDatatypeInfoToObject< GPUField<real_t> > ( bdId, numberOfGhostLayers );
      }

      IBlock * firstBlock =  & ( * bs->begin() );
      python_coupling::Dispatcher<GpuFields, Exporter_createMPIDatatypeInfoToObject > dispatcher( firstBlock );
      return dispatcher( bdId )( bdId, numberOfGhostLayers );
   }


   //===================================================================================================================
   //
   //  fieldCopy
   //
   //===================================================================================================================

   template<typename Field_T>
   void copyFieldToGpuDispatch(const shared_ptr<StructuredBlockStorage> & bs,
                               BlockDataID cpuFieldId, BlockDataID gpuFieldId, bool toGpu)
   {
      typedef cuda::GPUField<typename Field_T::value_type> GpuField;
      if(toGpu)
         cuda::fieldCpy<GpuField, Field_T>(bs, gpuFieldId, cpuFieldId);
      else
         cuda::fieldCpy<Field_T, GpuField>(bs, cpuFieldId, gpuFieldId);
   }
   FunctionExporterClass( copyFieldToGpuDispatch,
                          void( const shared_ptr<StructuredBlockStorage> &, BlockDataID, BlockDataID, bool ) );

   template< typename FieldTypes >
   void transferFields( const shared_ptr<StructuredBlockStorage> & bs,
                        const std::string & gpuFieldId, const std::string & cpuFieldId, bool toGpu)
   {
      if( bs->begin() == bs->end()) {
         return;
      };

      auto dstBdId = python_coupling::blockDataIDFromString( *bs, gpuFieldId );
      auto srcBdId = python_coupling::blockDataIDFromString( *bs, cpuFieldId );

      IBlock * firstBlock =  & ( * bs->begin() );
      python_coupling::Dispatcher<FieldTypes, Exporter_copyFieldToGpuDispatch> dispatcher( firstBlock );
      dispatcher( srcBdId )( bs, srcBdId, dstBdId, toGpu );
   }

   template< typename FieldTypes>
   void copyFieldToGpu(const shared_ptr<StructuredBlockStorage> & bs,
                       const std::string & gpuFieldId, const std::string & cpuFieldId)
   {
      transferFields<FieldTypes>(bs, gpuFieldId, cpuFieldId, true);
   }

   template< typename FieldTypes>
   void copyFieldToCpu(const shared_ptr<StructuredBlockStorage> & bs,
                       const std::string & gpuFieldId, const std::string & cpuFieldId)
   {
      transferFields<FieldTypes>(bs, gpuFieldId, cpuFieldId, false);
   }

} // namespace internal




template<typename GpuFields, typename CpuFields >
void exportModuleToPython()
{
   python_coupling::ModuleScope fieldModule( "cuda" );

   using namespace boost::python;

   python_coupling::for_each_noncopyable_type<GpuFields>( internal::GpuFieldExporter() );

   def( "createGpuField", &internal::createPythonGpuField<GpuFields>, ( ( arg("size")                     ),
                                                                         ( arg("type")                     ),
                                                                         ( arg("ghostLayers") = uint_t(1)  ),
                                                                         ( arg("layout")      = field::zyxf),
                                                                         ( arg("usePitchedMem") = true     )  ) );


   def( "addGpuFieldToStorage",  &internal::addToStorage<GpuFields>, ( ( arg("blocks")                  ),
                                                                        ( arg("name")                    ),
                                                                        ( arg("type")                    ),
                                                                        ( arg("fSize")       = 1         ),
                                                                        ( arg("ghostLayers") = uint_t(1) ),
                                                                        ( arg("layout")      = field::zyxf      ),
                                                                        ( arg("usePitchedMem") = object()  ) ) );

   def( "createMPIDatatypeInfo",&internal::createMPIDatatypeInfo<GpuFields>, ( arg("blocks"), arg("blockDataName"), arg("numberOfGhostLayers" ) =0 ) );
   def( "createPackInfo",       &internal::createPackInfo<GpuFields>,        ( arg("blocks"), arg("blockDataName"), arg("numberOfGhostLayers" ) =0 ) );

   def( "copyFieldToGpu", &internal::copyFieldToGpu<CpuFields>, (arg("blocks"), ("gpuFieldId"), ("cpuFieldId")));
   def( "copyFieldToCpu", &internal::copyFieldToCpu<CpuFields>, (arg("blocks"), ("gpuFieldId"), ("cpuFieldId")));
}





} // namespace cuda
} // namespace walberla


