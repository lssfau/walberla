//======================================================================================================================
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
//! \file WeightAssignmentFunctor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/amr/InfoCollection.h"

#include "blockforest/loadbalancing/PODPhantomData.h"

namespace walberla {
namespace pe {
namespace amr {

class WeightAssignmentFunctor
{
public:
   typedef walberla::blockforest::PODPhantomWeight<double>           PhantomBlockWeight;
   typedef walberla::blockforest::PODPhantomWeightPackUnpack<double> PhantomBlockWeightPackUnpackFunctor;

   ///Base weight due to allocated data structures. A weight of zero for blocks is dangerous as empty blocks might accumulate on one process!
   static const double baseWeight;

   WeightAssignmentFunctor( shared_ptr<InfoCollection>& ic ) : ic_(ic) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, boost::any > > & blockData, const PhantomBlockForest & )
   {
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         const PhantomBlock * block = it->first;
         //only change of one level is supported!
         WALBERLA_ASSERT_LESS( int_c(block->getLevel()) - int_c(block->getSourceLevel()), 2 );

         if (block->sourceBlockIsLarger())
         {
            auto infoIt = ic_->find( block->getId()/*.getFatherId()*/ );
            WALBERLA_ASSERT_UNEQUAL( infoIt, ic_->end() );
            it->second = PhantomBlockWeight( double_c(infoIt->second.numberOfLocalBodies) + baseWeight );
            continue;
         }

         if (block->sourceBlockHasTheSameSize())
         {
            auto infoIt = ic_->find( block->getId() );
            WALBERLA_ASSERT_UNEQUAL( infoIt, ic_->end() );
            it->second = PhantomBlockWeight( double_c(infoIt->second.numberOfLocalBodies) + baseWeight );
            continue;
         }

         if (block->sourceBlockIsSmaller())
         {
            double weight = 0;
            for (uint_t child = 0; child < 8; ++child)
            {
               blockforest::BlockID childId(block->getId(), child);
               auto childIt = ic_->find( childId );
               WALBERLA_ASSERT_UNEQUAL( childIt, ic_->end() );
               weight += double_c(childIt->second.numberOfLocalBodies);
            }
            it->second = PhantomBlockWeight( weight + baseWeight );
            continue;
         }
      }
   }

private:
   shared_ptr<InfoCollection> ic_;
};

}
}
}
