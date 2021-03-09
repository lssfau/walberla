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
//! \file DefaultBoundaryHandlingCollectionFactory.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "boundary/BoundaryHandlingCollection.h"

#include "lbm/boundary/factories/DefaultBoundaryHandling.h"
#include "lbm/boundary/factories/DefaultDiffusionBoundaryHandling.h"

#include <functional>

namespace walberla {
namespace lbm{

// TODO enable if depending on lattice model

template< typename LatticeModel_T, typename DiffusionLatticeModel_T, typename FlagField_T >
class DefaultBoundaryHandlingCollectionFactory
{
private:
   using DefaultBoundaryHandling_T = typename DefaultBoundaryHandlingFactory<LatticeModel_T, FlagField_T>::BoundaryHandling;
   using DefaultDiffusionBoundaryHandlingFactory_T = typename DefaultDiffusionBoundaryHandlingFactory<DiffusionLatticeModel_T, FlagField_T>::BoundaryHandling_T;

public:
   using BoundaryHandlingCollection_T = BoundaryHandlingCollection<FlagField_T, DefaultBoundaryHandling_T &, DefaultDiffusionBoundaryHandlingFactory_T &>;


private:
   static BoundaryHandlingCollection_T* createDefaultBoundaryHandlingCollectionFactory( 
      IBlock* const block, const StructuredBlockStorage* const /*bs*/, const BlockDataID& flagFieldID, const BlockDataID& handlingID, const BlockDataID& diffusionHandlingID )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      auto flagField = block->getData< FlagField_T >( flagFieldID );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      WALBERLA_ASSERT_NOT_NULLPTR( block->getData< DefaultBoundaryHandling_T                 >( handlingID          ) );
      WALBERLA_ASSERT_NOT_NULLPTR( block->getData< DefaultDiffusionBoundaryHandlingFactory_T >( diffusionHandlingID ) );

      DefaultBoundaryHandling_T                 & handling          = * block->getData< DefaultBoundaryHandling_T                 >( handlingID          );
      DefaultDiffusionBoundaryHandlingFactory_T & diffusionHandling = * block->getData< DefaultDiffusionBoundaryHandlingFactory_T >( diffusionHandlingID );
            
      return new BoundaryHandlingCollection_T( " Diffusion Boundary Handling Collection", flagField, handling, diffusionHandling );
   }

public:
   static BlockDataID addDefaultBoundaryHandlingCollectionToStorage(
      const shared_ptr< StructuredBlockStorage >& bs, const std::string & identifier, const BlockDataID& flagFieldID, const BlockDataID& handlingID, const BlockDataID& diffusionHandlingID )
   {
      auto func = std::bind( createDefaultBoundaryHandlingCollectionFactory, std::placeholders::_1, std::placeholders::_2, flagFieldID, handlingID, diffusionHandlingID );
      return bs->addStructuredBlockData< BoundaryHandlingCollection_T >( func, identifier );
   }

}; // class DefaultBoundaryHandlingCollectionFactory


} // namespace lbm
} // namespace walberla
