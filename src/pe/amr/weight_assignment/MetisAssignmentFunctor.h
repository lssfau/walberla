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
#include "domain_decomposition/PeriodicIntersectionVolume.h"

namespace walberla {
namespace pe {
namespace amr {

/**
 * Assignment functor for ParMetis based load balancing.
 *
 * Uses BlockInfo::computationalWeight as vertex weight.
 * Uses the overlap between neighboring domains as edge weight.
 */
class MetisAssignmentFunctor
{
public:
   ///Convenience typedef to be used as PhantomBlock weight in conjunction with this weight assignment functor.
   typedef blockforest::DynamicParMetisBlockInfo           PhantomBlockWeight;
   ///Convenience typdef for pack and unpack functions to be used in conjunction with PhantomBlockWeight.
   typedef blockforest::DynamicParMetisBlockInfoPackUnpack PhantomBlockWeightPackUnpackFunctor;

   MetisAssignmentFunctor( shared_ptr<InfoCollection>& ic, const real_t baseWeight = real_t(10.0) ) : ic_(ic), baseWeight_(baseWeight) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & forest )
   {
      const std::array< bool, 3 > periodic {{forest.getBlockForest().isPeriodic(0),
                  forest.getBlockForest().isPeriodic(1),
                  forest.getBlockForest().isPeriodic(2)}};
      const math::AABB domain     = forest.getBlockForest().getDomain();

      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         const PhantomBlock * block = it->first;
         //only change of one level is supported!
         WALBERLA_ASSERT_LESS( abs(int_c(block->getLevel()) - int_c(block->getSourceLevel())), 2 );

         //all information is provided by info collection
         auto infoIt         = ic_->find( block->getId() );
         WALBERLA_CHECK_UNEQUAL( infoIt, ic_->end() );
         const double weight = double_c( infoIt->second.computationalWeight ) + baseWeight_;
         blockforest::DynamicParMetisBlockInfo info( 0 );
         info.setVertexWeight( int64_c(weight) );
         info.setVertexSize( int64_c( weight ) );
         info.setVertexCoords( it->first->getAABB().center() );
         for( uint_t nb = uint_t(0); nb < it->first->getNeighborhoodSize(); ++nb )
         {
            const real_t dx(1.0);
            info.setEdgeWeight( it->first->getNeighborId(nb),
                                static_cast<blockforest::DynamicParMetisBlockInfo::weight_t>(
                                domain_decomposition::periodicIntersectionVolume( periodic,
                                                                                  domain,
                                                                                  it->first->getAABB(),
                                                                                  it->first->getNeighborAABB(nb).getExtended(dx))) );
         }
         it->second = info;
         continue;
      }
   }

   inline void   setBaseWeight( const double weight) { baseWeight_ = weight;}
   inline double getBaseWeight() const { return baseWeight_; }

private:
   shared_ptr< InfoCollection > ic_;

   ///Base weight due to allocated data structures. A weight of zero for blocks is dangerous as empty blocks might accumulate on one process!
   double baseWeight_ = 10.0;
};

}
}
}

