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
//! \file AddGPUFieldToStorage.impl.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "gpu/FieldCopy.h"

namespace walberla {
namespace gpu
{


   namespace internal
   {
      template< typename GPUField_T>
      GPUField_T * createGPUField( const IBlock * const block,
                                   const StructuredBlockStorage * const bs,
                                   uint_t ghostLayers,
                                   uint_t fSize,
                                   const field::Layout & layout,
                                   bool usePitchedMem )
      {
         return new GPUField_T( bs->getNumberOfXCells( *block ),
                                bs->getNumberOfYCells( *block ),
                                bs->getNumberOfZCells( *block ),
                                fSize, ghostLayers, layout, usePitchedMem );
      }

      template< typename Field_T>
      GPUField< typename Field_T::value_type> *
      createGPUFieldFromCPUField( const IBlock * const block,
                                  const StructuredBlockStorage * const,
                                  ConstBlockDataID cpuFieldID,
                                  bool usePitchedMem
                                )
      {
         using GPUField_T = GPUField<typename Field_T::value_type>;

         const Field_T * f = block->getData<Field_T>( cpuFieldID );
         auto gpuField = new GPUField_T( f->xSize(), f->ySize(), f->zSize(), f->fSize(),
                                         f->nrOfGhostLayers(), f->layout(), usePitchedMem );

         gpu::fieldCpy( *gpuField, *f );

         return gpuField;
      }

   }


   template< typename GPUField_T>
   BlockDataID addGPUFieldToStorage(const shared_ptr< StructuredBlockStorage >& bs,
                                    const std::string & identifier,
                                    uint_t fSize,
                                    const Layout layout,
                                    uint_t nrOfGhostLayers,
                                    bool usePitchedMem )
   {
      auto func = std::bind ( internal::createGPUField<GPUField_T>, std::placeholders::_1, std::placeholders::_2, nrOfGhostLayers, fSize, layout, usePitchedMem );
      return bs->addStructuredBlockData< GPUField_T >( func, identifier );
   }


   template< typename Field_T>
   BlockDataID addGPUFieldToStorage( const shared_ptr< StructuredBlockStorage > & bs,
                                     ConstBlockDataID cpuFieldID,
                                     const std::string & identifier,
                                     bool usePitchedMem )
   {
      auto func = [cpuFieldID, usePitchedMem](auto && PH1, auto && PH2) { return internal::createGPUFieldFromCPUField<Field_T>(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), cpuFieldID, usePitchedMem); };
      return bs->addStructuredBlockData< GPUField<typename Field_T::value_type> >( func, identifier );
   }



} // namespace gpu
} // namespace walberla
