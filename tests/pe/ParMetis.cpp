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
//! \file ParMetis.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include <blockforest/loadbalancing/DynamicParMetis.h>
#include <core/debug/TestSubsystem.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Sample.h>

using namespace walberla;

class ReGrid
{
public:
   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & /*forest*/ )
   {
      std::for_each( minTargetLevels.begin(),
                     minTargetLevels.end(),
                     [](auto& pair){pair.second = pair.first->getLevel() + 1;} );
   }

};

class MetisAssignmentFunctor
{
public:

   using PhantomBlockWeight = blockforest::DynamicParMetisBlockInfo;
   using PhantomBlockWeightPackUnpackFunctor = blockforest::DynamicParMetisBlockInfoPackUnpack;

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
   {
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         const auto&   corner = it->first->getAABB().maxCorner();
         const int weight     = int_c( corner[0] + corner[1] + corner[2] );
         blockforest::DynamicParMetisBlockInfo info(0);
         info.setVertexWeight( weight );
         info.setVertexSize( weight );
         info.setVertexCoords( it->first->getAABB().center() );
         for( uint_t nb = uint_t(0); nb < it->first->getNeighborhoodSize(); ++nb )
         {
            info.setEdgeWeight(it->first->getNeighborId(nb), int64_c(weight) );
         }
         it->second = info;
      }
   }
};

int parmetisTest(const std::string& algorithm,
                 const std::string& weightsToUse,
                 const std::string& edgeSource)
{
   walberla::MPIManager::instance()->resetMPI();
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT("****** " << algorithm << " | " << weightsToUse << " | " << edgeSource);

   // create forest
   auto forest = blockforest::createBlockForest( math::AABB(0,0,0,40,40,40),
                                                 Vector3<uint_t>(4, 4, 4),
                                                 Vector3<bool>(false, false, false),
                                                 64,
                                                 0);
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }

   forest->setRefreshMinTargetLevelDeterminationFunction( ReGrid() );

   auto assFunc = MetisAssignmentFunctor();
   forest->setRefreshPhantomBlockDataAssignmentFunction( assFunc );
   forest->setRefreshPhantomBlockDataPackFunction( MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
   forest->setRefreshPhantomBlockDataUnpackFunction( MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

   auto alg     = blockforest::DynamicParMetis::stringToAlgorithm(    algorithm );
   auto vWeight = blockforest::DynamicParMetis::stringToWeightsToUse( weightsToUse );
   auto eWeight = blockforest::DynamicParMetis::stringToEdgeSource(   edgeSource );

   auto prepFunc = blockforest::DynamicParMetis( alg, vWeight, eWeight );
   forest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );

   forest->refresh();

   math::Sample numBlocks;
   numBlocks.castToRealAndInsert( forest->size() );
   numBlocks.mpiGatherRoot();
   WALBERLA_LOG_INFO_ON_ROOT("#blocks: " << numBlocks.format() );

   int weight = 0;
   for (const auto& block : *forest)
   {
      const auto&   corner = block.getAABB().maxCorner();
      weight    += int_c( corner[0] + corner[1] + corner[2] );
   }
   math::Sample weightSample;
   weightSample.castToRealAndInsert( weight );
   weightSample.mpiGatherRoot();
   WALBERLA_LOG_INFO_ON_ROOT("#weights: " << weightSample.format() );

   return EXIT_SUCCESS;
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::vector<std::string> algs = {"PART_GEOM", "PART_GEOM_KWAY", "PART_KWAY", "ADAPTIVE_REPART", "REFINE_KWAY"};
   std::vector<std::string> wtu  = {"NO_WEIGHTS", "EDGE_WEIGHTS", "VERTEX_WEIGHTS", "BOTH_WEIGHTS"};
   std::vector<std::string> es   = {"EDGES_FROM_FOREST", "EDGES_FROM_EDGE_WEIGHTS"};

   for (const auto& a : algs)
      for (const auto& w : wtu)
         for (const auto& e : es)
         {
            parmetisTest(a,w,e);
         }

   return EXIT_SUCCESS;
}
