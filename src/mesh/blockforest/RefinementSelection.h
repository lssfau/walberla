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
//! \file RefinementSelection.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh/MeshOperations.h"

#include "blockforest/SetupBlockForest.h"
#include "blockforest/Types.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/Optional.h"

#include <random>
#include <vector>

namespace walberla {
namespace mesh {

template< typename DistanceObject >
class RefinementSelection
{
public:
   typedef typename DistanceObject::Scalar Scalar;

   RefinementSelection( const shared_ptr< DistanceObject > & distanceObject, const uint_t level, const real_t distance, real_t maxError ) : distanceObject_( distanceObject ), level_( level ), distance_( distance ), maxError_( maxError )  {  }

   void operator()( blockforest::SetupBlockForest & forest ) const;

private:

   shared_ptr< DistanceObject > distanceObject_;
   
   uint_t level_;
   real_t distance_;
   real_t maxError_;
};

template< typename DistanceObject >
inline RefinementSelection< DistanceObject > makeRefinementSelection( const shared_ptr< DistanceObject > & distanceObject, const uint_t level, const real_t distance, real_t maxError )
{
   return RefinementSelection< DistanceObject >( distanceObject, level, distance, maxError );
}

template< typename DistanceObject >
inline void walberla::mesh::RefinementSelection<DistanceObject>::operator()( blockforest::SetupBlockForest & forest ) const
{
   std::vector< blockforest::SetupBlock* > blocks;
   forest.getBlocks( blocks );

   const uint_t numBlocks    = uint_c( blocks.size() );
   const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );
   const uint_t chunkSize    = uint_c( std::ceil( real_c( numBlocks ) / real_c( numProcesses ) ) );

   const uint_t rank    = uint_c( MPIManager::instance()->rank() );
   const int chunkBegin = int_c( rank * chunkSize );
   const int chunkEnd   = std::min( int_c( ( rank + 1 ) * chunkSize ), int_c( numBlocks ) );

   std::vector<size_t> shuffle( numBlocks );
   for( size_t i = 0; i < shuffle.size(); ++i )
   {
      shuffle[i] = i;
   }

   //old compilers might need this
   //std::srand( 42 );
   //std::random_shuffle( shuffle.begin(), shuffle.end() );

   std::mt19937 g( 42 );
   std::shuffle( shuffle.begin(), shuffle.end(), g );

   std::vector<uint8_t> refine( numBlocks, uint_t( 0 ) );
   
   #ifdef _OPENMP
   #pragma omp parallel for schedule( dynamic )
   #endif
   for( int i = chunkBegin; i < chunkEnd; ++i )
   {
      const size_t ii = numeric_cast<size_t>( i );
      const uint_t blockLevel = blocks[ shuffle[ii] ]->getLevel();
      if( blockLevel >= level_ )
         continue;

      walberla::optional< bool > intersects = intersectsSurface( *distanceObject_, blocks[ shuffle[ii] ]->getAABB(), maxError_, distance_ );
      if( !intersects || intersects.value() )
         refine[ shuffle[ ii ] ] = uint8_t( 1 );
   }

   allReduceInplace( refine, mpi::LOGICAL_OR );

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      blocks[i]->setMarker( refine[ i ] ==  uint8_t( 1 ) );
   }
}

} // namespace mesh
} // namespace walberla


