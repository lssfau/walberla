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
//! \file MetisAssignmentFunctor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/amr/InfoCollection.h"

#include "blockforest/loadbalancing/DynamicParMetis.h"

namespace walberla {
namespace pe {
namespace amr {

class MetisAssignmentFunctor
{
public:

   typedef blockforest::DynamicParMetisBlockInfo           PhantomBlockWeight;
   typedef blockforest::DynamicParMetisBlockInfoPackUnpack PhantomBlockWeightPackUnpackFunctor;

   MetisAssignmentFunctor( shared_ptr<InfoCollection>& ic, const real_t baseWeight = real_t(10.0) ) : ic_(ic), baseWeight_(baseWeight) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
   {
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         const double weight     = double_c( ic_->find( it->first->getId() )->second.numberOfLocalBodies ) + baseWeight_;
         //const double commWeight = double_c( edgeWeightFactor * (double_c(ic_->find( it->first->getId() )->second.numberOfShadowBodies) + baseWeight_)) + 1;
         blockforest::DynamicParMetisBlockInfo info( 0 );
         info.setVertexWeight( int64_c(weight) );
         info.setVertexSize( int64_c( weight ) );
         info.setVertexCoords( it->first->getAABB().center() );
         for( uint_t nb = uint_t(0); nb < it->first->getNeighborhoodSize(); ++nb )
         {
            info.setEdgeWeight(it->first->getNeighborId(nb), int64_c(edgeWeight_) );
         }
         it->second = info;
      }
   }

   inline void   setBaseWeight( const double weight) { baseWeight_ = weight;}
   inline double getBaseWeight() const { return baseWeight_; }

   inline void   setEdgeWeight( const double weight) { edgeWeight_ = weight;}
   inline double getEdgeWeight() const { return edgeWeight_; }

private:
   shared_ptr< InfoCollection > ic_;

   ///Base weight due to allocated data structures. A weight of zero for blocks is dangerous as empty blocks might accumulate on one process!
   double baseWeight_ = 10.0;
   double edgeWeight_ = 1.0;
};

}
}
}

