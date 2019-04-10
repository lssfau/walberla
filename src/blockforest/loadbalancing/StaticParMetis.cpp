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
//! \file StaticParMetis.cpp
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "StaticParMetis.h"

#include "core/load_balancing/ParMetisWrapper.h"

#include "core/logging/Logging.h"

#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Gather.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/Tokenizing.h"

#include "core/timing/Timer.h"

#include "core/StringUtility.h"

namespace walberla {
namespace blockforest {

template< typename T >
T * ptr( std::vector<T> & v )
{
   if( v.empty() )
      return nullptr;
   else
      return &( v.front() );
}

using idx_t = uint_t;


uint_t StaticLevelwiseParMetis::operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t /*perProcessMemoryLimit*/ ) const
{
   WcTimer globalTimer, parmetisTimer;

   int numRunnerProcesses = MPIManager::instance()->numProcesses(); // Number of processes running ParMetis (!= number of processes the partition is computed for)
   int rank = MPIManager::instance()->rank();

   for( uint_t level = forest.getMinLevel(); level <= forest.getMaxLevel(); ++level )
   {
      globalTimer.start();
      WALBERLA_LOG_INFO_ON_ROOT( "Running static ParMetis partitioning on blocks of level " << level );

      std::vector< SetupBlock* > blocks;
      forest.getBlocks( blocks, level );

      const uint_t numBlocks = blocks.size();

      if( numBlocks < numberOfProcesses )
      {
         for( uint_t i = 0; i < numBlocks; ++i )
            blocks[i]->assignTargetProcess( i );

         WALBERLA_LOG_INFO_ON_ROOT( "Less blocks on level " << level << " (" << numBlocks << ") than processes (" << numberOfProcesses << "). Using simple load balancing."  );

         globalTimer.end();
         continue;
      }

      for( uint_t i = 0; i < uint_c( blocks.size() ); ++i )
      {
         blocks[i]->setIndex( i );
      }

      const uint_t chunkSize  = uint_c( std::ceil( real_c( numBlocks ) / real_c( numRunnerProcesses ) ) );
      const uint_t chunkBegin = std::min( chunkSize * uint_c( rank ), numBlocks );
      const uint_t chunkEnd   = std::min( chunkSize * uint_c( rank + 1 ), numBlocks );

      std::vector<int64_t> vtxdist;
      vtxdist.reserve( uint_c(numRunnerProcesses) + uint_t(1) );
      for( uint_t i = 0; i < uint_c(numRunnerProcesses); ++i )
         vtxdist.push_back( int64_c( std::min( i * chunkSize, numBlocks ) ) );
      vtxdist.push_back( int64_t( forest.getNumberOfBlocks( level ) ) );

      std::vector<int64_t> adjncy, xadj, vwgt, adjwgt;
      std::vector<double> xyz;
      std::vector< BlockPair > blockPairs;

      for( uint_t i = chunkBegin; i < chunkEnd; ++i )
      {
         const SetupBlock & block = *blocks[i];

         xadj.push_back( int64_c( adjncy.size() ) );


         for( auto nit = block.getNeighborhood().begin(); nit != block.getNeighborhood().end(); ++nit )
         {
            if( (*nit)->getLevel() != level )
               continue; // ignore neighbor blocks on other levels

            adjncy.push_back( int64_c( (*nit)->getIndex() ) );

            if(weightsToUse_ == PARMETIS_EDGE_WEIGHTS || weightsToUse_ == PARMETIS_BOTH_WEIGHTS)
            {
               blockPairs.emplace_back( blocks[i], *nit );
            }
         }

         vwgt.push_back( int64_c( block.getWorkload() ) );
         Vector3<real_t> center = block.getAABB().center();
         xyz.push_back( center[0] );
         xyz.push_back( center[1] );
         xyz.push_back( center[2] );
      }
      xadj.push_back( int64_c( adjncy.size() ) );

      if( weightsToUse_ == PARMETIS_EDGE_WEIGHTS || weightsToUse_ == PARMETIS_BOTH_WEIGHTS )
      {
         WALBERLA_ASSERT_EQUAL( adjncy.size(), blockPairs.size() );
         adjwgt.resize( blockPairs.size(), int64_t(1) );
         commWeightFunction_( blockPairs, adjwgt );

         if( adjwgt.empty() )
            adjwgt.push_back( int64_t(0) ); // dummy value to circumvent dumb NULL pointer check of ParMetis
      }

      WALBERLA_ASSERT_EQUAL( vtxdist.size(), uint_c(numRunnerProcesses) + uint_t( 1 ) );
      WALBERLA_ASSERT_EQUAL( uint_c(vtxdist[uint_c(rank)]), chunkBegin );
      WALBERLA_ASSERT_EQUAL( uint_c(vtxdist[uint_c(rank + 1)]), chunkEnd );
      WALBERLA_ASSERT_EQUAL( xadj.size(), (chunkEnd - chunkBegin) + 1 );
      WALBERLA_ASSERT_EQUAL( vwgt.size(), chunkEnd - chunkBegin );

      std::vector<int64_t> part( chunkEnd - chunkBegin, int64_c( MPIManager::instance()->rank() ) );

      int64_t wgtflag = weightsToUse_;
      int64_t numflag = 0; // C-style ordering
      int64_t ncon = 1; // Number of constraints
      int64_t ndims = 3; // Number of dimensions
      double ubvec[] = { real_t( 1.05 ) }; // imbalance tolerance
      int64_t nparts = int64_c( numberOfProcesses ); // number of subdomains
      MPI_Comm comm = MPIManager::instance()->comm();
      std::vector<double> tpwgts( uint_c(nparts * ncon), 1.0 / double_c( nparts ) ); // vertex weight fraction that is stored in a subdomain
      int64_t options[] = { int64_t( 1 ), int64_t( 0 ), int64_t( 23 ), int64_t( 1 ) };

      // add dummy element to circumvent null pointer check if less blocks than processes
      adjncy.resize( std::max( adjncy.size(), size_t(1) ) );
      part.resize( std::max( part.size(), size_t(1) ) );
      vwgt.resize( std::max( vwgt.size(), size_t(1) ) );

      int64_t edgecut = 0;
      int metisResult = core::METIS_ERROR;

      switch( algorithm_ )
      {
      case PARMETIS_PART_GEOM_KWAY:
         parmetisTimer.start();
         metisResult = core::ParMETIS_V3_PartGeomKway( ptr( vtxdist ), ptr( xadj ), ptr( adjncy ), ptr( vwgt ), ptr( adjwgt ), &wgtflag, &numflag, &ndims, ptr( xyz ), &ncon, &nparts, ptr( tpwgts ), ubvec, options, &edgecut, ptr( part ), &comm );
         parmetisTimer.end();
         break;
      case PARMETIS_PART_KWAY:
         parmetisTimer.start();
         metisResult = core::ParMETIS_V3_PartKway( ptr( vtxdist ), ptr( xadj ), ptr( adjncy ), ptr( vwgt ), ptr( adjwgt ), &wgtflag, &numflag, &ncon, &nparts, ptr( tpwgts ), ubvec, options, &edgecut, ptr( part ), &comm );
         parmetisTimer.end();
         break;
      }

      if( metisResult != core::METIS_OK )
         WALBERLA_ABORT("ParMetis failed!");

      std::vector< int64_t > parts = mpi::allGatherv( part, comm );

      WALBERLA_ASSERT_EQUAL( parts.size(), numBlocks );

      for( uint_t i = 0; i < numBlocks; ++i )
         blocks[i]->assignTargetProcess( uint_c( parts[i] ) );

      globalTimer.end();
      WALBERLA_LOG_INFO_ON_ROOT( "ParMetis partitioning finished for level " << level << " successfully after " << globalTimer.last() << " s (ParMetis took " << parmetisTimer.last() << " s = " << std::setprecision(2) << parmetisTimer.last() / globalTimer.last() * 100.0 << "%) with an edge cut of " << edgecut );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "ParMetis partitioning finished after " << globalTimer.total() << " s (ParMetis took " << parmetisTimer.total() << " s = " << std::setprecision(2) << parmetisTimer.total() / globalTimer.total() * 100.0 << "%)" );

   //count number of used processes
   std::vector<bool> processUsed( numberOfProcesses, false );
   for(auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt)
   {
      processUsed[ blockIt->getTargetProcess() ] = true;
   }

   return uint_c(std::count( processUsed.begin(), processUsed.end(), true ));
}

StaticLevelwiseParMetis::Algorithm StaticLevelwiseParMetis::stringToAlgorithm( std::string s )
{
   string_to_upper( s );
   string_trim( s );

   if( s == "PART_GEOM_KWAY" )
      return PARMETIS_PART_GEOM_KWAY;
   else if( s == "PART_KWAY" )
      return PARMETIS_PART_KWAY;
   else
      WALBERLA_ABORT( "Illegal ParMetis algorithm specified! Valid choices are: \"PART_GEOM_KWAY\" or \"PART_KWAY\"." );
}


StaticLevelwiseParMetis::WeightsToUse StaticLevelwiseParMetis::stringToWeightsToUse( std::string s )
{
   string_to_upper( s );
   string_trim( s );

   if( s == "NO_WEIGHTS" )
      return PARMETIS_NO_WEIGHTS;
   else if( s == "EDGE_WEIGHTS" )
      return PARMETIS_EDGE_WEIGHTS;
   else if( s == "VERTEX_WEIGHTS" )
      return PARMETIS_VERTEX_WEIGHTS;
   else if( s == "BOTH_WEIGHTS" )
      return PARMETIS_BOTH_WEIGHTS;
   else
      WALBERLA_ABORT( "Illegal ParMetis weights usage specified! Valid choices are: \"NO_WEIGHTS\", \"EDGE_WEIGHTS\", \"VERTEX_WEIGHTS\", or \"BOTH_WEIGHTS\"." );
}


} // namespace blockforest
} // namespace walberla
