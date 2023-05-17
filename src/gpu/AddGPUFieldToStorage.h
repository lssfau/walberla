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
//! \file AddGPUFieldToStorage.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "GPUField.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace gpu
{



   //*******************************************************************************************************************
   /*! Adds a gpu::GPUField to a StructuredBlockStorage
   *
   *  - Similar to walberla::field::addToStorage() functions
   *  - created field is uninitialized
   */
   //*******************************************************************************************************************
   template< typename GPUField_T>
   BlockDataID addGPUFieldToStorage(const shared_ptr< StructuredBlockStorage >& bs,
                                    const std::string & identifier,
                                    uint_t fSize,
                                    const Layout layout = fzyx,
                                    uint_t nrOfGhostLayers = 1,
                                    bool usePitchedMem = true );



   //*******************************************************************************************************************
   /*! Adds a gpu::GPUField to a StructuredBlockStorage using data from a CPU field
   *
   *  - adds a GPU field to a StructuredBlockStorage using a CPU field
   *  - sizes, number of ghostlayers and layout are the same as the CPU field
   *  - GPU field is initialized with the data currently stored in the CPU field
   *  @tparam Field_T  type of the CPU field, the created GPUField will be of type gpu::GPUField<Field_T::value_type>
   */
   //*******************************************************************************************************************
   template< typename Field_T>
   BlockDataID addGPUFieldToStorage( const shared_ptr< StructuredBlockStorage > & bs,
                                     ConstBlockDataID cpuFieldID,
                                     const std::string & identifier,
                                     bool usePitchedMem = true );



} // namespace gpu
} // namespace walberla


#include "AddGPUFieldToStorage.impl.h"
