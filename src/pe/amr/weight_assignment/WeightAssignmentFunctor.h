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

#include "blockforest/loadbalancing/InfoCollection.h"
#include "blockforest/loadbalancing/PODPhantomData.h"
#include "blockforest/PhantomBlockForest.h"

namespace walberla {
namespace pe {
namespace amr {

/**
 * General assignment functor for load balancing.
 *
 * Uses BlockInfo::computationalWeight as PhantomBlockWeight
 */
class WeightAssignmentFunctor
{
public:
   ///Convenience typedef to be used as PhantomBlock weight in conjunction with this weight assignment functor.
   using PhantomBlockWeight = walberla::blockforest::PODPhantomWeight<double>;
   ///Convenience typedef for pack and unpack functions to be used in conjunction with PhantomBlockWeight.
   using PhantomBlockWeightPackUnpackFunctor = walberla::blockforest::PODPhantomWeightPackUnpack<double>;

   WeightAssignmentFunctor( shared_ptr<blockforest::InfoCollection>& ic, const real_t baseWeight = real_t(10.0) ) : ic_(ic), baseWeight_(baseWeight) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
   {
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         const PhantomBlock * block = it->first;
         //only change of one level is supported!
         WALBERLA_CHECK_LESS( abs(int_c(block->getLevel()) - int_c(block->getSourceLevel())), 2 );

         auto infoIt = ic_->find( block->getId()/*.getFatherId()*/ );
         WALBERLA_CHECK_UNEQUAL( infoIt, ic_->end() );
         it->second = PhantomBlockWeight( double_c(infoIt->second.computationalWeight) + baseWeight_ );
      }
   }

   inline void   setBaseWeight( const double weight) { baseWeight_ = weight;}
   inline double getBaseWeight() const { return baseWeight_; }

private:
   shared_ptr<blockforest::InfoCollection> ic_;

   ///Base weight due to allocated data structures. A weight of zero for blocks is dangerous as empty blocks might accumulate on one process!
   double baseWeight_ = real_t(10.0);
};

}
}
}
