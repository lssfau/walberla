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
//! \file DynamicParMetis.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/PhantomBlockForest.h"

#include "core/debug/Debug.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIWrapper.h"

#include <cmath>
#include <map>

namespace walberla {
namespace blockforest {

std::pair<uint_t, uint_t> getBlockSequenceRange( const PhantomBlockForest & phantomForest, MPI_Comm comm );

std::map< blockforest::BlockID, uint_t > getBlockIdToSequenceMapping( const PhantomBlockForest & phantomForest, const std::pair<uint_t, uint_t> & blockSequenceRange, MPI_Comm comm );


class DynamicParMetis
{
public:
   enum Algorithm    { PARMETIS_PART_GEOM_KWAY, PARMETIS_PART_KWAY, PARMETIS_ADAPTIVE_REPART, PARMETIS_REFINE_KWAY };
   enum WeightsToUse { PARMETIS_NO_WEIGHTS = 0, PARMETIS_EDGE_WEIGHTS = 1, PARMETIS_VERTEX_WEIGHTS = 2, PARMETIS_BOTH_WEIGHTS = 3 };
   enum EdgeSource   { PARMETIS_EDGES_FROM_FOREST, PARMETIS_EDGES_FROM_EDGE_WEIGHTS };

   DynamicParMetis( const Algorithm algorithm = PARMETIS_PART_GEOM_KWAY,
                    const WeightsToUse weightsToUse = PARMETIS_BOTH_WEIGHTS,
                    const EdgeSource edgeSource = PARMETIS_EDGES_FROM_EDGE_WEIGHTS )
      : algorithm_( algorithm ), weightsToUse_( weightsToUse ), edgeSource_( edgeSource ) {}

   bool operator()( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                    std::set< uint_t > & processesToRecvFrom,
                    const PhantomBlockForest & phantomForest,
                    const uint_t iteration ) const;


   void   setipc2redist(double val) {ipc2redist_ = val;}
   double getipc2redist() const {return ipc2redist_;}

   bool edgeWeightsUsed()   const { return ( weightsToUse_ == PARMETIS_EDGE_WEIGHTS   ) || ( weightsToUse_ == PARMETIS_BOTH_WEIGHTS ); }
   bool vertexWeightsUsed() const { return ( weightsToUse_ == PARMETIS_VERTEX_WEIGHTS ) || ( weightsToUse_ == PARMETIS_BOTH_WEIGHTS ); }
   bool vertexSizeUsed()    const { return algorithm_ == PARMETIS_ADAPTIVE_REPART; }

   static Algorithm    stringToAlgorithm( std::string s );
   static WeightsToUse stringToWeightsToUse( std::string s );
   static EdgeSource   stringToEdgeSource( std::string s );

   std::string algorithmToString() const;
   std::string weightsToUseToString() const;
   std::string edgeSourceToString() const;

protected:
   Algorithm algorithm_;
   WeightsToUse weightsToUse_;
   EdgeSource edgeSource_;

   double ipc2redist_ = real_t( 1000.0 ); ///< compute repartitioning with low edge cut (set lower (down to 0.000001) to get minimal repartitioning )
};

class DynamicParMetisBlockInfo
{
public:

   typedef int64_t weight_t;
   typedef int64_t vsize_t;

   DynamicParMetisBlockInfo( const weight_t vertexWeight )
      : vertexWeight_(vertexWeight), vertexSize_(1)
   { }

   DynamicParMetisBlockInfo( mpi::RecvBuffer & buffer )
   {
      buffer >> vertexWeight_ >> vertexSize_ >> vertexCoords_ >> edgeWeights_;
   }

   weight_t getVertexWeight() const { return vertexWeight_; }
   void     setVertexWeight( const weight_t vertexWeight ) { vertexWeight_ = vertexWeight; }

   vsize_t getVertexSize() const { return vertexSize_; }
   void setVertexSize( const vsize_t size ) { vertexSize_ = size; }

   const Vector3<double> & getVertexCoords() const { WALBERLA_ASSERT( !std::isnan(vertexCoords_[0]) && !std::isnan(vertexCoords_[1]) && !std::isnan(vertexCoords_[2]) ); return vertexCoords_; }
   void setVertexCoords( const Vector3<double> & p ) { vertexCoords_ = p; }

   void setEdgeWeight( const blockforest::BlockID & blockId, const weight_t edgeWeight ) { edgeWeights_[blockId] = edgeWeight; }
   void setEdgeWeights( const std::map< blockforest::BlockID, weight_t > & edgeWeights ) { edgeWeights_ = edgeWeights; }

   const std::map< blockforest::BlockID, weight_t > & getEdgeWeights() const { return edgeWeights_; }

   void toBuffer( mpi::SendBuffer & buffer )
   {
      buffer << vertexWeight_ << vertexSize_ << vertexCoords_ << edgeWeights_;
   }

private:

   weight_t vertexWeight_; /// Defines the weight of a block
   vsize_t vertexSize_;    /// Defines the cost of rebalancing a block
   Vector3<double> vertexCoords_ = Vector3<double>(std::numeric_limits<double>::signaling_NaN()); /// Defines where the block is located in space. Needed by some ParMetis algorithms
   std::map< blockforest::BlockID, weight_t > edgeWeights_; /// Defines the cost of communication with other blocks
};


struct DynamicParMetisBlockInfoPackUnpack
{
   void operator()( mpi::SendBuffer & buffer, const PhantomBlock & block )
   {
      block.getData< DynamicParMetisBlockInfo >().toBuffer( buffer );
   }

   void operator()( mpi::RecvBuffer & buffer, const PhantomBlock &, walberla::any & data )
   {
      data = DynamicParMetisBlockInfo( buffer );
   }
};







} // namespace blockforest
} // namespace walberla
