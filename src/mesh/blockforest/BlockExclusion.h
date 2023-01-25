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
//! \file BlockExclusion.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh_common/MeshOperations.h"

#include "blockforest/SetupBlockForest.h"

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include <random>
#include <vector>

namespace walberla::mesh {

template< typename DistanceObject >
class ExcludeMeshExterior
{
public:
   ExcludeMeshExterior( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) : distanceObject_( distanceObject ), maxError_( maxError ) { }

   void operator()( std::vector<uint8_t> & excludeBlock, const blockforest::SetupBlockForest::RootBlockAABB & aabb ) const;

private:
   shared_ptr< DistanceObject > distanceObject_;
   real_t maxError_;
};

template< typename DistanceObject >
class ExcludeMeshInteriorRefinement
{
 public:
   ExcludeMeshInteriorRefinement( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) : distanceObject_( distanceObject ), maxError_( maxError ) { }

   bool operator()( const blockforest::SetupBlock & block ) const;

 private:
   shared_ptr< DistanceObject > distanceObject_;
   real_t maxError_;
};


template< typename DistanceObject >
class ExcludeMeshInterior
{
public:
   ExcludeMeshInterior( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) : distanceObject_( distanceObject ), maxError_( maxError ) { }

   void operator()( std::vector<uint8_t> & excludeBlock, const blockforest::SetupBlockForest::RootBlockAABB & aabb ) const;

private:
   shared_ptr< DistanceObject > distanceObject_;
   real_t maxError_;
};


template< typename DistanceObject >
ExcludeMeshExterior<DistanceObject> makeExcludeMeshExterior( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) { return ExcludeMeshExterior<DistanceObject>( distanceObject, maxError ); }


template< typename DistanceObject >
ExcludeMeshInterior<DistanceObject> makeExcludeMeshInterior( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) { return ExcludeMeshInterior<DistanceObject>( distanceObject, maxError ); }


template< typename DistanceObject >
ExcludeMeshInteriorRefinement<DistanceObject> makeExcludeMeshInteriorRefinement( const shared_ptr< DistanceObject > & distanceObject, const real_t maxError ) { return ExcludeMeshInteriorRefinement<DistanceObject>( distanceObject, maxError ); }

template< typename DistanceObject >
void walberla::mesh::ExcludeMeshExterior<DistanceObject>::operator()( std::vector<uint8_t> & excludeBlock, const blockforest::SetupBlockForest::RootBlockAABB & aabb ) const
{
      const uint_t numBlocks    = uint_c( excludeBlock.size() );
      const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );
      const uint_t chunkSize    = uint_c( std::ceil( real_c( numBlocks ) / real_c( numProcesses ) ) );

      const uint_t rank       = uint_c( MPIManager::instance()->rank() );
      const size_t chunkBegin = rank * chunkSize;
      const size_t chunkEnd   = std::min( ( rank + 1 ) * chunkSize, numBlocks );

      std::vector<size_t> shuffle( excludeBlock.size() );
      for( size_t i = 0; i < excludeBlock.size(); ++i )
      {
         shuffle[i] = i;
      }

      //old compilers might need this
      //std::srand( 42 );
      //std::random_shuffle( shuffle.begin(), shuffle.end() );

      std::mt19937 g( 42 );
      std::shuffle( shuffle.begin(), shuffle.end(), g );

      std::fill( excludeBlock.begin(), excludeBlock.end(), uint8_t( 0 ) );

      #ifdef _OPENMP
      #pragma omp parallel for schedule( dynamic )
      #endif
      for( size_t i = chunkBegin; i < chunkEnd; ++i )
      {
         auto intersectionDefined = isIntersecting( *distanceObject_, aabb( shuffle[i] ), maxError_ );
         if( intersectionDefined && !intersectionDefined.value() )
            excludeBlock[ shuffle[i] ] = uint8_t( 1 );
      }

      allReduceInplace( excludeBlock, mpi::LOGICAL_OR );
}


template< typename DistanceObject >
void walberla::mesh::ExcludeMeshInterior<DistanceObject>::operator()( std::vector<uint8_t> & excludeBlock, const blockforest::SetupBlockForest::RootBlockAABB & aabb ) const
{
   const uint_t numBlocks    = uint_c( excludeBlock.size() );
   const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );
   const uint_t chunkSize    = uint_c( std::ceil( real_c( numBlocks ) / real_c( numProcesses ) ) );

   const uint_t rank    = uint_c( MPIManager::instance()->rank() );
   const int chunkBegin = int_c( rank * chunkSize );
   const int chunkEnd   = std::min( int_c( ( rank + 1 ) * chunkSize ), int_c( numBlocks ) );

   std::vector<size_t> shuffle( excludeBlock.size() );
   for( size_t i = 0; i < excludeBlock.size(); ++i )
   {
      shuffle[i] = i;
   }

   //old compilers might need this
   //std::srand( 42 );
   //std::random_shuffle( shuffle.begin(), shuffle.end() );

   std::mt19937 g( 42 );
   std::shuffle( shuffle.begin(), shuffle.end(), g );

   std::fill( excludeBlock.begin(), excludeBlock.end(), uint8_t( 0 ) );

#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic )
#endif
   for( int i = chunkBegin; i < chunkEnd; ++i )
   {
      auto is = numeric_cast<size_t>( i );
      auto fullCoveringAABBDefined = fullyCoversAABB( *distanceObject_, aabb( shuffle[is] ), maxError_ );
      if( fullCoveringAABBDefined && fullCoveringAABBDefined.value() )
         excludeBlock[ shuffle[is] ] = uint8_t( 1 );
   }

   allReduceInplace( excludeBlock, mpi::LOGICAL_OR );
}

template< typename DistanceObject >
bool walberla::mesh::ExcludeMeshInteriorRefinement<DistanceObject>::operator()( const blockforest::SetupBlock & block ) const
{
   const AABB aabb = block.getAABB();
   auto fullCoveringAABBDefined = fullyCoversAABB( *distanceObject_, aabb, maxError_ );

   return static_cast<bool>(fullCoveringAABBDefined && fullCoveringAABBDefined.value());
}

} // namespace walberla::mesh

