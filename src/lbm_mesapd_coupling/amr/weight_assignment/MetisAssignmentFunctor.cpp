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
//! \file MetisAssignmentFunctor.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "MetisAssignmentFunctor.h"

namespace walberla {
namespace lbm_mesapd_coupling {
namespace amr {

void MetisAssignmentFunctor::operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData,
      const PhantomBlockForest & /*phantomBlockForest*/)
{
   for (auto &it : blockData) {
      const PhantomBlock * block = it.first;
      //only change of one level is supported!
      WALBERLA_ASSERT_LESS( std::abs(int_c(block->getLevel()) - int_c(block->getSourceLevel())), 2 );

      std::vector<int64_t> metisVertexWeights(weightEvaluationFctVector_.size());

      for( auto con = uint_t(0); con < weightEvaluationFctVector_.size(); ++con )
      {
         real_t vertexWeight = std::max(weightEvaluationFctVector_[con](block), blockBaseWeight_);

         int64_t metisVertexWeight = int64_c( weightMultiplicator_ * vertexWeight );

         WALBERLA_ASSERT_GREATER(metisVertexWeight, int64_t(0));
         metisVertexWeights[con] = metisVertexWeight;
      }

      blockforest::DynamicParMetisBlockInfo info( metisVertexWeights );

      info.setVertexCoords(it.first->getAABB().center() );

      real_t blockVolume = it.first->getAABB().volume();
      real_t approximateEdgeLength = std::cbrt( blockVolume );

      int64_t faceNeighborWeight = int64_c(approximateEdgeLength * approximateEdgeLength ); //common face
      int64_t edgeNeighborWeight = int64_c( approximateEdgeLength ); //common edge
      int64_t cornerNeighborWeight = int64_c( 1 ); //common corner

      int64_t vertexSize = int64_c(blockVolume);
      WALBERLA_ASSERT_GREATER(vertexSize, int64_t(0));
      info.setVertexSize( vertexSize );

      for( const uint_t idx : blockforest::getFaceNeighborhoodSectionIndices() )
      {
         for( auto nb = uint_t(0); nb < it.first->getNeighborhoodSectionSize(idx); ++nb )
         {
            auto neighborBlockID = it.first->getNeighborId(idx,nb);
            info.setEdgeWeight(neighborBlockID, faceNeighborWeight );
         }
      }

      for( const uint_t idx : blockforest::getEdgeNeighborhoodSectionIndices() )
      {
         for( auto nb = uint_t(0); nb < it.first->getNeighborhoodSectionSize(idx); ++nb )
         {
            auto neighborBlockID = it.first->getNeighborId(idx,nb);
            info.setEdgeWeight(neighborBlockID, edgeNeighborWeight );
         }
      }

      for( const uint_t idx : blockforest::getCornerNeighborhoodSectionIndices() )
      {
         for( auto nb = uint_t(0); nb < it.first->getNeighborhoodSectionSize(idx); ++nb )
         {
            auto neighborBlockID = it.first->getNeighborId(idx,nb);
            info.setEdgeWeight(neighborBlockID, cornerNeighborWeight );
         }
      }

      it.second = info;

   }
}


} // namespace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla
