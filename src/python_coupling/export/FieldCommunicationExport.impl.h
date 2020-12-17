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
//! \file CommunicationExport.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#   include "blockforest/StructuredBlockForest.h"
#   include "python_coupling/helper/BlockStorageExportHelpers.h"

#   include "field/communication/PackInfo.h"
#   include "field/communication/StencilRestrictedPackInfo.h"
#   include "field/communication/UniformMPIDatatypeInfo.h"

#   include "python_coupling/helper/MplHelpers.h"

#   include "stencil/D2Q9.h"
#   include "stencil/D3Q15.h"
#   include "stencil/D3Q19.h"
#   include "stencil/D3Q27.h"
#   include "stencil/D3Q7.h"

#  include <typeinfo>

#  include "pybind11/pybind11.h"

namespace walberla
{
namespace field
{
namespace internal
{
namespace py = pybind11;

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

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> )
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;
      typedef GhostLayerField<T, F_SIZE> GlField_T;
      IBlock * firstBlock =  & ( * blocks_->begin() );
      if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
      {
         if ( numberOfGhostLayers_ > 0  )
         {
            resultPackInfo_ = py::cast(make_shared< field::communication::PackInfo< GlField_T > >(fieldId_, numberOfGhostLayers_));
         }
         else
         {
            resultPackInfo_ = py::cast(make_shared< field::communication::PackInfo< GlField_T > >(fieldId_));
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


template<typename... FieldTypes>
static py::object PackInfoWrapper(const shared_ptr<StructuredBlockForest> & blocks,
                                  const std::string & name, uint_t numberOfGhostLayers )
{
   BlockDataID fieldID = python_coupling::blockDataIDFromString( *blocks, name );

   if ( blocks->begin() == blocks->end() ) {
      // if no blocks are on this field an arbitrary PackInfo can be returned
      return py::cast( make_shared< field::communication::PackInfo<GhostLayerField<real_t,1>> >( fieldID, numberOfGhostLayers ) );
   }

   PackInfoExporter exporter(blocks, fieldID, numberOfGhostLayers);
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( std::ref(exporter) );
   if ( ! exporter.getResultPackInfo() ) {
      throw py::value_error("Failed to create PackInfo");
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

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> )
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;
      typedef GhostLayerField<T, F_SIZE> GlField_T;
      IBlock * firstBlock =  & ( * blocks_->begin() );
      if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
      {
         if ( numberOfGhostLayers_ > 0  )
            resultMPIDatatypeInfo_ =  py::cast( make_shared< field::communication::UniformMPIDatatypeInfo<GlField_T> >( fieldId_, numberOfGhostLayers_ ) );
         else
            resultMPIDatatypeInfo_ =  py::cast( make_shared< field::communication::UniformMPIDatatypeInfo<GlField_T> >( fieldId_ ) );

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


template<typename... FieldTypes>
static py::object UniformMPIDatatypeInfoWrapper(const shared_ptr<StructuredBlockForest> & blocks,
                                  const std::string & name, uint_t numberOfGhostLayers )
{
   BlockDataID fieldID = python_coupling::blockDataIDFromString( *blocks, name );

   if ( blocks->begin() == blocks->end() ) {
      // if no blocks are on this field an arbitrary PackInfo can be returned
      return py::cast( make_shared< field::communication::UniformMPIDatatypeInfo<GhostLayerField<real_t,1>> >( fieldID, numberOfGhostLayers ) );
   }

   UniformMPIDatatypeInfoExporter exporter(blocks, fieldID, numberOfGhostLayers);
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( std::ref(exporter) );
   if ( ! exporter.getResultUniformMPIDatatype() ) {
      throw py::value_error("Failed to create UniformMPIDatatype");
   }
   else {
      return exporter.getResultUniformMPIDatatype();
   }
}

//===================================================================================================================
//
//  exportStencilRestrictedPackInfo
//
//===================================================================================================================

template< typename T >
void exportStencilRestrictedPackInfo(py::module_& m)
{
   using field::communication::StencilRestrictedPackInfo;
   {
      typedef StencilRestrictedPackInfo< GhostLayerField< T, 9 >, stencil::D2Q9 > Pi;
      py::class_< Pi, shared_ptr< Pi >, walberla::communication::UniformPackInfo >(m, "StencilRestrictedPackInfo_D2Q9");
   }
   {
      typedef StencilRestrictedPackInfo< GhostLayerField< T, 7 >, stencil::D3Q7 > Pi;
      py::class_< Pi, shared_ptr< Pi >, walberla::communication::UniformPackInfo >(m, "StencilRestrictedPackInfo_D3Q7");
   }
   {
      typedef StencilRestrictedPackInfo< GhostLayerField< T, 15 >, stencil::D3Q15 > Pi;
      py::class_< Pi, shared_ptr< Pi >, walberla::communication::UniformPackInfo >(m,
                                                                                   "StencilRestrictedPackInfo_D3Q15");
   }
   {
      typedef StencilRestrictedPackInfo< GhostLayerField< T, 19 >, stencil::D3Q19 > Pi;
      py::class_< Pi, shared_ptr< Pi >, walberla::communication::UniformPackInfo >(m,
                                                                                   "StencilRestrictedPackInfo_D3Q19");
   }
   {
      typedef StencilRestrictedPackInfo< GhostLayerField< T, 27 >, stencil::D3Q27 > Pi;
      py::class_< Pi, shared_ptr< Pi >, walberla::communication::UniformPackInfo >(m,
                                                                                   "StencilRestrictedPackInfo_D3Q27");
   }
}

} // namespace internal

namespace py = pybind11;
using namespace pybind11::literals;

template< typename... FieldTypes >
void exportCommunicationClasses(py::module_& m)
{
   py::module_ m2 = m.def_submodule("field", "Field Extension of the waLBerla python bindings");
   internal::exportStencilRestrictedPackInfo< real_t >(m2);

   m2.def(
      "createPackInfo",
      [](const shared_ptr<StructuredBlockForest> & blocks,
         const std::string & blockDataName, uint_t numberOfGhostLayers ) {
        return internal::PackInfoWrapper< FieldTypes... >(blocks, blockDataName, numberOfGhostLayers);
      },
      "blocks"_a, "blockDataName"_a, "numberOfGhostLayers"_a = uint_t(0));

   m2.def(
      "createMPIDatatypeInfo",
      [](const shared_ptr<StructuredBlockForest> & blocks,
         const std::string & blockDataName, uint_t numberOfGhostLayers ) {
        return internal::UniformMPIDatatypeInfoWrapper< FieldTypes... >(blocks, blockDataName, numberOfGhostLayers);
      },
      "blocks"_a, "blockDataName"_a, "numberOfGhostLayers"_a = uint_t(0));
}

} // namespace field
} // namespace walberla

#endif // WALBERLA_BUILD_WITH_PYTHON
