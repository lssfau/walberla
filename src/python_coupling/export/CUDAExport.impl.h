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
//! \file CUDAExport.impl.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "core/logging/Logging.h"

#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/FieldCopy.h"
#include "cuda/GPUField.h"
#include "cuda/communication/GPUPackInfo.h"

#include "field/AddToStorage.h"
#include "field/communication/UniformMPIDatatypeInfo.h"

#include "python_coupling/PythonWrapper.h"
#include "python_coupling/helper/MplHelpers.h"

namespace walberla {
namespace cuda {



namespace internal {
using namespace pybind11::literals;
   //===================================================================================================================
   //
   //  Field export
   //
   //===================================================================================================================

   template<typename GpuField_T>
   uint64_t gpufield_ptr(const GpuField_T & gpuField)
   {
      return reinterpret_cast<uint64_t>(gpuField.pitchedPtr().ptr);
      // return gpuField.pitchedPtr();
   }

   template<typename GpuField_T>
   std::string gpufield_dtypeStr(const GpuField_T & )
   {
      return std::string(field::internal::PythonFormatString<typename GpuField_T::value_type>::get());
   }

   struct GpuFieldExporter
   {
      GpuFieldExporter(py::module_& m) : m_(m) {}
      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> ) const
      {

         typedef typename GpuField_T::value_type T;
         std::string data_type_name = field::internal::PythonFormatString<T>::get();

         std::string class_name = "GpuField_" + data_type_name;
         py::class_<GpuField_T, shared_ptr<GpuField_T>>(m_, class_name.c_str() )
            .def_property_readonly("layout",              &field::internal::field_layout            < GpuField_T > )
            .def_property_readonly("size",                &field::internal::field_size              < GpuField_T > )
            .def_property_readonly("sizeWithGhostLayers", &field::internal::field_sizeWithGhostLayer< GpuField_T > )
            .def_property_readonly("allocSize",           &field::internal::field_allocSize         < GpuField_T > )
            .def_property_readonly("strides",             &field::internal::field_strides           < GpuField_T > )
            .def_property_readonly("offsets",             &field::internal::field_offsets           < GpuField_T > )
            .def_property_readonly("ptr",                 &gpufield_ptr                             < GpuField_T > )
            .def_property_readonly("dtypeStr",            &gpufield_dtypeStr                        < GpuField_T > )
            .def_property_readonly("isPitchedMem",        &GpuField_T::isPitchedMem )
            .def("swapDataPointers",    &field::internal::field_swapDataPointers  < GpuField_T > )
            .def_property_readonly("nrOfGhostLayers",     &GpuField_T::nrOfGhostLayers )
            .def("cloneUninitialized",  &GpuField_T::cloneUninitialized, py::return_value_policy::copy)
            ;


         using field::communication::PackInfo;
         using communication::GPUPackInfo;
         std::string GpuFieldPackInfoName = "GpuFieldPackInfo_" + data_type_name;
         py::class_< GPUPackInfo<GpuField_T>, shared_ptr< GPUPackInfo<GpuField_T> >, walberla::communication::UniformPackInfo>(m_, GpuFieldPackInfoName.c_str() );

         using field::communication::UniformMPIDatatypeInfo;
         std::string GpuFieldMPIDataTypeInfoName = "GpuFieldMPIDataTypeInfo_" + data_type_name;
         py::class_< UniformMPIDatatypeInfo<GpuField_T>, shared_ptr< UniformMPIDatatypeInfo<GpuField_T> >, walberla::communication::UniformMPIDatatypeInfo>(m_, GpuFieldMPIDataTypeInfoName.c_str() );

      }
      const py::module_& m_;
   };


   //===================================================================================================================
   //
   //  addToStorage
   //
   //===================================================================================================================

   class AddToStorageExporter
   {
   public:
      AddToStorageExporter( const shared_ptr<StructuredBlockForest> & blocks, const std::string & name,
                            py::object &dtype, uint_t fs, uint_t gl, Layout layout,
                            bool usePitchedMem )
         : blocks_( blocks ), name_( name ), dtype_(dtype), fs_( fs ),
           gl_(gl),layout_( layout), usePitchedMem_(usePitchedMem), found_(true)
      {}

      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         typedef typename GpuField_T::value_type T;

         if(python_coupling::isCppEqualToPythonType<T>(py::cast<std::string>(dtype_.attr("__name__"))))
         {
            addGPUFieldToStorage< GPUField< T > >(blocks_, name_, fs_, layout_, gl_, usePitchedMem_);
         }
      }

      bool successful() const { return found_; }
   private:
      shared_ptr< StructuredBlockForest > blocks_;
      std::string name_;
      py::object dtype_;
      uint_t fs_;
      uint_t gl_;
      Layout layout_;
      bool usePitchedMem_;
      bool found_;
   };

   template<typename... GpuFields>
   void addToStorage( const shared_ptr<StructuredBlockForest> & blocks, const std::string & name, py::object &dtype,
                      uint_t fs, uint_t gl, Layout layout, bool usePitchedMem )
   {
      namespace py = pybind11;
      auto result = make_shared<py::object>();
      AddToStorageExporter exporter( blocks, name, dtype, fs, gl, layout, usePitchedMem );
      python_coupling::for_each_noncopyable_type<GpuFields...>( std::ref(exporter) );
   }


   //===================================================================================================================
   //
   //  createPackInfo Export
   //
   //===================================================================================================================

   class PackInfoExporter
   {
    public:
      PackInfoExporter(const shared_ptr<StructuredBlockForest> & blocks, BlockDataID fieldId, uint_t numberOfGhostLayers)
         : blocks_(blocks), fieldId_(fieldId), numberOfGhostLayers_( numberOfGhostLayers )
      {}

      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         using cuda::communication::GPUPackInfo;

         IBlock * firstBlock =  & ( * blocks_->begin() );
         if( firstBlock->isDataClassOrSubclassOf<GpuField_T>(fieldId_) )
         {
            if ( numberOfGhostLayers_ > 0  )
            {
               resultPackInfo_ = py::cast(make_shared< GPUPackInfo< GpuField_T > >(fieldId_, numberOfGhostLayers_));
            }
            else
            {
               resultPackInfo_ = py::cast(make_shared< GPUPackInfo< GpuField_T > >(fieldId_));
            }
         }
      }
      py::object getResultPackInfo()
      {
         return resultPackInfo_;
      }

    private:
      py::object resultPackInfo_;
      shared_ptr< StructuredBlockStorage > blocks_;
      BlockDataID fieldId_;
      uint_t numberOfGhostLayers_;
   };


   template<typename... GpuField_T>
   static py::object PackInfoWrapper(const shared_ptr<StructuredBlockForest> & blocks,
                                     const std::string & name, uint_t numberOfGhostLayers )
   {
      using cuda::communication::GPUPackInfo;
      BlockDataID fieldID = python_coupling::blockDataIDFromString( *blocks, name );

      if ( blocks->begin() == blocks->end() ) {
         // if no blocks are on this field an arbitrary PackInfo can be returned
         return py::cast( make_shared< GPUPackInfo<GPUField<int8_t>> >( fieldID, numberOfGhostLayers ) );
      }

      PackInfoExporter exporter(blocks, fieldID, numberOfGhostLayers);
      python_coupling::for_each_noncopyable_type< GpuField_T... >  ( std::ref(exporter) );
      if ( ! exporter.getResultPackInfo() ) {
         throw py::value_error("Failed to create GPU PackInfo");
      }
      else {
         return exporter.getResultPackInfo();
      }
   }

   //===================================================================================================================
   //
   //  createMPIDatatypeInfo
   //
   //===================================================================================================================

   class UniformMPIDatatypeInfoExporter
   {
    public:
      UniformMPIDatatypeInfoExporter(const shared_ptr<StructuredBlockForest> & blocks, BlockDataID fieldId, uint_t numberOfGhostLayers)
         : blocks_(blocks), fieldId_(fieldId), numberOfGhostLayers_( numberOfGhostLayers )
      {}

      template< typename GpuField_T>
      void operator() ( python_coupling::NonCopyableWrap<GpuField_T> )
      {
         using field::communication::UniformMPIDatatypeInfo;
         IBlock * firstBlock =  & ( * blocks_->begin() );
         if( firstBlock->isDataClassOrSubclassOf<GpuField_T>(fieldId_) )
         {
            if ( numberOfGhostLayers_ > 0  )
               resultMPIDatatypeInfo_ =  py::cast( make_shared< UniformMPIDatatypeInfo<GpuField_T> >( fieldId_, numberOfGhostLayers_ ) );
            else
               resultMPIDatatypeInfo_ =  py::cast( make_shared< UniformMPIDatatypeInfo<GpuField_T> >( fieldId_ ) );

         }
      }
      py::object getResultUniformMPIDatatype()
      {
         return resultMPIDatatypeInfo_;
      }

    private:
      py::object resultMPIDatatypeInfo_;
      shared_ptr< StructuredBlockStorage > blocks_;
      BlockDataID fieldId_;
      uint_t numberOfGhostLayers_;
   };


   template<typename... GpuField_T>
   static py::object UniformMPIDatatypeInfoWrapper(const shared_ptr<StructuredBlockForest> & blocks,
                                                   const std::string & name, uint_t numberOfGhostLayers )
   {
      using field::communication::UniformMPIDatatypeInfo;
      BlockDataID fieldID = python_coupling::blockDataIDFromString( *blocks, name );

      if ( blocks->begin() == blocks->end() ) {
         // if no blocks are on this field an arbitrary PackInfo can be returned
         return py::cast( make_shared< UniformMPIDatatypeInfo<GPUField<int8_t>> >( fieldID, numberOfGhostLayers ) );
      }

      UniformMPIDatatypeInfoExporter exporter(blocks, fieldID, numberOfGhostLayers);
      python_coupling::for_each_noncopyable_type< GpuField_T... >  ( std::ref(exporter) );
      if ( ! exporter.getResultUniformMPIDatatype() ) {
         throw py::value_error("Failed to create GPU UniformMPIDatatype");
      }
      else {
         return exporter.getResultUniformMPIDatatype();
      }
   }

   //===================================================================================================================
   //
   //  fieldCopy
   //
   //===================================================================================================================

class copyFieldToGpuDispatchExporter
{
 public:
   copyFieldToGpuDispatchExporter( const shared_ptr<StructuredBlockForest> & blocks,
                                   BlockDataID gpuFieldId, BlockDataID cpuFieldId, bool toGPU)
      : blocks_( blocks ), gpuFieldId_( gpuFieldId ), cpuFieldId_(cpuFieldId), toGPU_( toGPU )
   {}

   template< typename CpuField_T>
   void operator() ( python_coupling::NonCopyableWrap<CpuField_T> )
   {
      typedef cuda::GPUField<typename CpuField_T::value_type> GpuField_T;
      IBlock * firstBlock =  & ( * blocks_->begin() );

      if(firstBlock->isDataClassOrSubclassOf< CpuField_T > ( cpuFieldId_ ) )
      {
         if(toGPU_)
           cuda::fieldCpy<GpuField_T, CpuField_T>(blocks_, gpuFieldId_, cpuFieldId_);
         else
           cuda::fieldCpy<CpuField_T, GpuField_T>(blocks_, cpuFieldId_, gpuFieldId_);
      }
   }
 private:
   shared_ptr< StructuredBlockForest > blocks_;
   BlockDataID gpuFieldId_;
   BlockDataID cpuFieldId_;
   bool toGPU_;
};

template<typename... CpuFields>
void copyFieldToGPU(const shared_ptr< StructuredBlockForest > & blocks, const std::string & gpuFieldName,
                    const std::string & cpuFieldName, bool toGPU )
{
   namespace py = pybind11;
   auto result = make_shared<py::object>();

   BlockDataID gpuFieldId = python_coupling::blockDataIDFromString( *blocks, gpuFieldName );
   BlockDataID cpuFieldId = python_coupling::blockDataIDFromString( *blocks, cpuFieldName );

   copyFieldToGpuDispatchExporter exporter( blocks, gpuFieldId, cpuFieldId, toGPU );
   python_coupling::for_each_noncopyable_type<CpuFields...>( std::ref(exporter) );
}
} // namespace internal


using namespace pybind11::literals;

template<typename... GpuFields>
void exportModuleToPython(py::module_ &m)
{
   py::module_ m2 = m.def_submodule("cuda", "Cuda Extension of the waLBerla python bindings");

   python_coupling::for_each_noncopyable_type<GpuFields...>( internal::GpuFieldExporter(m2) );

   m2.def(
      "addGpuFieldToStorage",
      [](const shared_ptr< StructuredBlockForest > & blocks, const std::string & name, py::object &dtype, uint_t fSize,
         bool usePitchedMem, uint_t ghostLayers, Layout layout) {
        return internal::addToStorage<GpuFields...>(blocks, name, dtype, fSize, ghostLayers, layout, usePitchedMem);
      },
      "blocks"_a, "name"_a, "dtype"_a, "fSize"_a=1, "usePitchedMem"_a=false, "ghostLayers"_a=uint(1), "layout"_a=fzyx);

   m2.def(
      "createPackInfo",
      [](const shared_ptr<StructuredBlockForest> & blocks,
         const std::string & blockDataName, uint_t numberOfGhostLayers ) {
         return internal::PackInfoWrapper< GpuFields... >(blocks, blockDataName, numberOfGhostLayers);
      },
      "blocks"_a, "blockDataName"_a, "numberOfGhostLayers"_a = uint_t(0));

   m2.def(
      "createMPIDatatypeInfo",
      [](const shared_ptr<StructuredBlockForest> & blocks,
         const std::string & blockDataName, uint_t numberOfGhostLayers ) {
        return internal::UniformMPIDatatypeInfoWrapper< GpuFields... >(blocks, blockDataName, numberOfGhostLayers);
      },
      "blocks"_a, "blockDataName"_a, "numberOfGhostLayers"_a = uint_t(0));

}

template<typename... CpuFields >
void exportCopyFunctionsToPython(py::module_ &m)
{
     py::module_ m2 = m.def_submodule("cuda", "Cuda Extension of the waLBerla python bindings");

   m2.def(
      "copyFieldToGpu",
      [](const shared_ptr< StructuredBlockForest > & blocks, const std::string & gpuFieldName, const std::string & cpuFieldName) {
        return internal::copyFieldToGPU<CpuFields...>(blocks, gpuFieldName, cpuFieldName, true);
      },
      "blocks"_a, "gpuFieldName"_a, "cpuFieldName"_a);

   m2.def(
      "copyFieldToCpu",
      [](const shared_ptr< StructuredBlockForest > & blocks, const std::string & gpuFieldName, const std::string & cpuFieldName) {
        return internal::copyFieldToGPU<CpuFields...>(blocks, gpuFieldName, cpuFieldName, false);
      },
      "blocks"_a, "gpuFieldName"_a, "cpuFieldName"_a);
}




} // namespace cuda
} // namespace walberla


