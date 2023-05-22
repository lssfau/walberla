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
//! \file AddToStorage.h
//! \ingroup lbm_generated
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/uid/SUID.h"

#include "gpu/GPUWrapper.h"
#include "gpu/FieldCopy.h"

#include "field/blockforest/BlockDataHandling.h"

#include "GPUPdfField.h"

namespace walberla::lbm_generated
{

namespace internal
{

template< typename LatticeStorageSpecification_T>
GPUPdfField< LatticeStorageSpecification_T > * createGPUPdfField( const IBlock * const block,
                           const StructuredBlockStorage * const bs,
                           const LatticeStorageSpecification_T& storageSpecification,
                           const uint_t ghostLayers,
                           const field::Layout & layout,
                           const bool usePitchedMem )
{
   using GPUField_T = GPUPdfField< LatticeStorageSpecification_T >;

   auto gpuField = new GPUField_T(bs->getNumberOfXCells( *block ),
                                  bs->getNumberOfYCells( *block ),
                                  bs->getNumberOfZCells( *block ),
                                  storageSpecification, ghostLayers,
                                  layout, usePitchedMem);

   return gpuField;
}

template< typename Field_T, typename LatticeStorageSpecification_T >
GPUPdfField< LatticeStorageSpecification_T >*
   createGPUPdfFieldFromCPUPdfField(const IBlock* const block, const StructuredBlockStorage* const,
                                    const LatticeStorageSpecification_T& storageSpecification,
                                    ConstBlockDataID cpuFieldID, const bool usePitchedMem, const bool copyCPUField = true)
{
   using GPUField_T = GPUPdfField< LatticeStorageSpecification_T >;

   const Field_T* f = block->getData< Field_T >(cpuFieldID);

   auto gpuField = new GPUField_T(f->xSize(), f->ySize(), f->zSize(), storageSpecification, f->nrOfGhostLayers(),
                                  f->layout(), usePitchedMem);

   if (copyCPUField)
      gpu::fieldCpy(*gpuField, *f);

   return gpuField;
}

} // namespace internal

template< typename GPUField_T, typename LatticeStorageSpecification_T >
BlockDataID addGPUPdfFieldToStorage(const shared_ptr< StructuredBlockStorage >& bs,
                                    const std::string & identifier,
                                    const LatticeStorageSpecification_T& storageSpecification,
                                    const Layout layout = fzyx,
                                    const uint_t nrOfGhostLayers = 1,
                                    const bool usePitchedMem = true )
{

   auto func = std::bind(internal::createGPUPdfField< LatticeStorageSpecification_T >,
                         std::placeholders::_1, std::placeholders::_2, storageSpecification, nrOfGhostLayers, layout, usePitchedMem);
   return bs->addStructuredBlockData< GPUPdfField< LatticeStorageSpecification_T > >(func, identifier);
}

template< typename Field_T, typename LatticeStorageSpecification_T >
BlockDataID addGPUPdfFieldToStorage(const shared_ptr< StructuredBlockStorage >& bs, ConstBlockDataID cpuFieldID,
                                    const LatticeStorageSpecification_T& storageSpecification,
                                    const std::string& identifier, const bool usePitchedMem = true, const bool copyCPUField = true)
{
   auto func = std::bind(internal::createGPUPdfFieldFromCPUPdfField< Field_T, LatticeStorageSpecification_T >,
                         std::placeholders::_1, std::placeholders::_2, storageSpecification, cpuFieldID, usePitchedMem, copyCPUField);
   return bs->addStructuredBlockData< GPUPdfField< LatticeStorageSpecification_T > >(func, identifier);
}

} // namespace walberla::lbm_generated