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
//! \file DynamicDiffusive.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "NoPhantomData.h"
#include "blockforest/BlockForest.h"
#include "blockforest/PhantomBlockForest.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/math/Random.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include <map>
#include <set>
#include <type_traits>
#include <vector>



namespace walberla {
namespace blockforest {


/**
 *  This class implements a diffusive algorithm for load balancing.
 *
 *  All algorithms are implemented to work levelwise. Load balancing with levels ignored is possible
 *  by specifying levelwise = false in the constructor.
**/
template< typename PhantomData_T >
class DynamicDiffusionBalance
{
public:

   enum Mode { DIFFUSION_PUSH, DIFFUSION_PULL, DIFFUSION_PUSHPULL };

   DynamicDiffusionBalance( const uint_t maxIterations, const uint_t flowIterations, const bool levelwise = true ) :
      mode_( DIFFUSION_PUSHPULL ), maxIterations_( maxIterations ),
      defineProcessWeightLimitByMultipleOfMaxBlockWeight_( true ), checkForEarlyAbort_( true ), abortThreshold_( 1.0 ),
      adaptOutflowWithGlobalInformation_( true ), adaptInflowWithGlobalInformation_( true ),
      flowIterations_( flowIterations ), flowIterationsIncreaseStart_( maxIterations ), flowIterationsIncrease_( 0.0 ),
      regardConnectivity_( true ), disregardConnectivityStart_( maxIterations ), outflowExceedFactor_( 1.0 ), inflowExceedFactor_( 1.0 ),
      levelwise_(levelwise)
   {}
   
   void setMode( const Mode mode ) { mode_ = mode; }
   Mode getMode() const { return mode_; }
   
   void setMaxIterations( const uint_t maxIterations ) { maxIterations_ = maxIterations; }
   uint_t getMaxIterations() const { return maxIterations_; }
   
   void defineProcessWeightLimitByMultipleOfMaxBlockWeight( const bool b ) { defineProcessWeightLimitByMultipleOfMaxBlockWeight_ = b; }
   bool defineProcessWeightLimitByMultipleOfMaxBlockWeight() const { return defineProcessWeightLimitByMultipleOfMaxBlockWeight_; }

   void checkForEarlyAbort( const bool check, const double abortThreshold = 1.0 )
   {
      checkForEarlyAbort_ = check;
      abortThreshold_ = abortThreshold;
   }
   bool checkForEarlyAbort() const { return checkForEarlyAbort_; }
   
   void adaptOutflowWithGlobalInformation( const bool adapt ) { adaptOutflowWithGlobalInformation_ = adapt; }
   bool adaptOutflowWithGlobalInformation() const { return adaptOutflowWithGlobalInformation_; }
   
   void adaptInflowWithGlobalInformation( const bool adapt ) { adaptInflowWithGlobalInformation_ = adapt; }
   bool adaptInflowWithGlobalInformation() const { return adaptInflowWithGlobalInformation_; }
   
   void setFlowIterations( const uint_t flowIterations ) { flowIterations_ = flowIterations; }
   uint_t getFlowIterations() const { return flowIterations_; }
   
   void setDynamicFlowIterationsIncrease( const uint_t startingBalanceIteration, const double increasePerBalanceIteration )
   {
      flowIterationsIncreaseStart_ = startingBalanceIteration;
      flowIterationsIncrease_ = increasePerBalanceIteration;
   }
   
   void regardConnectivity( const bool c ) { regardConnectivity_ = c; }
   bool regardConnectivity() const { return regardConnectivity_; }
   
   void disregardConnectivity( const uint_t startingBalanceIteration ) { disregardConnectivityStart_ = startingBalanceIteration; }
   
   void setOutflowExceedFactor( const double f ) { outflowExceedFactor_ = f; }
   double getOutflowExceedFactor() const { return outflowExceedFactor_; }
   
   void setInflowExceedFactor( const double f ) { inflowExceedFactor_ = f; }
   double getInflowExceedFactor() const { return inflowExceedFactor_; }
   
   bool operator()( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                    std::set< uint_t > & processesToRecvFrom,
                    const PhantomBlockForest & phantomForest,
                    const uint_t iteration );

private:

   double weight( const PhantomBlock * block ) const
   {
      return std::is_same< PhantomData_T, NoPhantomData >::value ? 1.0 :
               numeric_cast< double >( block->template getData< PhantomData_T >().weight() );
   }

   Mode mode_;
   
   uint_t maxIterations_;
   
   bool defineProcessWeightLimitByMultipleOfMaxBlockWeight_; // only evaluated when checkForEarlyAbort_ == true or adaptOutflowWithGlobalInformation_ == true
   bool checkForEarlyAbort_;
   double abortThreshold_; // only evaluated when checkForEarlyAbort_ == true
   bool adaptOutflowWithGlobalInformation_;
   bool adaptInflowWithGlobalInformation_;
   
   uint_t flowIterations_;
   uint_t flowIterationsIncreaseStart_;
   double flowIterationsIncrease_;
   
   bool regardConnectivity_;
   uint_t disregardConnectivityStart_;
   double outflowExceedFactor_;
   double inflowExceedFactor_;
   
   math::IntRandom< uint_t > random_;

   /// All gets for levels are wrapped like
   /// \code
   /// levelwise_ ? getCorrectLevel() : 0
   /// \endcode
   ///
   /// This allows to use the same algorithm for levelwise balancing as well as for balancing without levels.
   bool levelwise_;
};



template< typename PhantomData_T >
bool DynamicDiffusionBalance< PhantomData_T >::operator()( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                    std::set< uint_t > & processesToRecvFrom,
                                                                    const PhantomBlockForest & phantomForest,
                                                                    const uint_t iteration )
{
   auto & blockforest = phantomForest.getBlockForest();
   auto & neighborhood = phantomForest.getNeighborhood();
   
   const uint_t levels = levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1);
   
   // determine process weight (for every level)

   std::vector< double > processWeight( levels, 0.0 );
   std::vector< double > avgProcessWeight( levels, 0.0 );
   std::vector< double > maxBlockWeight( levels, 0.0 );
   std::vector< double > processWeightLimit( levels, 0.0 );
   
   //fill processWeight with total weight per level
   //find maxBlockWeight per level
   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
   {
      const uint_t level = levelwise_ ? it->first->getLevel() : uint_t(0);
      const auto blockWeight = weight( it->first );
      WALBERLA_ASSERT_LESS( level, levels );
      WALBERLA_CHECK_GREATER_EQUAL( blockWeight, 0.0 );
      processWeight[ level ] += blockWeight;
      maxBlockWeight[ level ] = std::max( blockWeight, maxBlockWeight[ level ] );
   }
   
   // determine avg. process weight and max. block weight (for every level)
   // and use this data to check if an early abort is possible (because the load is already balanced)
   
   if( checkForEarlyAbort_ || adaptOutflowWithGlobalInformation_ )
   {
      avgProcessWeight = processWeight;
      mpi::allReduceInplace( avgProcessWeight, mpi::SUM, MPIManager::instance()->comm() );
      const double numProcesses = double_c( MPIManager::instance()->numProcesses() );
      for( auto it = avgProcessWeight.begin(); it != avgProcessWeight.end(); ++it )
         *it /= numProcesses;
      mpi::allReduceInplace( maxBlockWeight, mpi::MAX, MPIManager::instance()->comm() );
   }
   
   for( uint_t l = uint_t(0); l != levels; ++l )
   {
      if( defineProcessWeightLimitByMultipleOfMaxBlockWeight_ ) // calculate process weight limit as multiple of max block weight
      {
         processWeightLimit[l] = maxBlockWeight[l];
         while( processWeightLimit[l] < avgProcessWeight[l] )
            processWeightLimit[l] += maxBlockWeight[l];
      }
      else
      {
         processWeightLimit[l] = std::max( avgProcessWeight[l], maxBlockWeight[l] );
      }
   }
   
   uint_t levelsToProcess( levels );
   std::vector< bool > processLevel( levels, true ); // determines whether the corresponding level must be processed
   if( checkForEarlyAbort_ )
   {
      for( uint_t l = uint_t(0); l != levels; ++l )
         processLevel[l] = ( processWeight[l] > processWeightLimit[l] * abortThreshold_ );
      mpi::allReduceInplace( processLevel, mpi::BITWISE_OR, MPIManager::instance()->comm() );
      
      // If no level needs processing, we can perform an early abort of the load balancing step.
      
      bool abort( true );
      for( uint_t l = uint_t(0); l != levels; ++l )
      {
         abort = ( abort && !(processLevel[l]) );
         if( !(processLevel[l]) )
            levelsToProcess -= uint_t(1);
      }
      if( abort )
         return false;
   }
   WALBERLA_ASSERT_GREATER( levelsToProcess, uint_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( levelsToProcess, levels );

   // alpha exchange

   std::map< mpi::MPIRank, mpi::MPISize > alphaRanksToRecvFrom;
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      alphaRanksToRecvFrom[ static_cast< mpi::MPIRank >(*n) ] = mpi::BufferSizeTrait<double>::size;

   mpi::BufferSystem alphaBufferSystem( MPIManager::instance()->comm(), 1708 ); // dynamicdiffusion = 100 121 110 097 109 105 099 100 105 102 102 117 115 105 111 110
   alphaBufferSystem.setReceiverInfo( alphaRanksToRecvFrom );

   std::map< uint_t, double > alpha; //process rank -> alpha
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
   {
      WALBERLA_ASSERT( alpha.find(*n) == alpha.end() );
      alpha[*n] = 1.0 / ( double_c( neighborhood.size() ) + 1.0 );
   }

   for( auto rank = alphaRanksToRecvFrom.begin(); rank != alphaRanksToRecvFrom.end(); ++rank )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c(rank->first), blockforest.getProcess() );
      WALBERLA_ASSERT( alpha.find( uint_c(rank->first) ) != alpha.end() );
      alphaBufferSystem.sendBuffer( rank->first ) << alpha[ uint_c(rank->first) ];
   }

   alphaBufferSystem.sendAll();

   for( auto recvIt = alphaBufferSystem.begin(); recvIt != alphaBufferSystem.end(); ++recvIt )
   {
      uint_t np = uint_c( recvIt.rank() );
      WALBERLA_ASSERT_UNEQUAL( np, blockforest.getProcess() );
      WALBERLA_ASSERT( alpha.find( np ) != alpha.end() );
      double a( 0.0 );
      recvIt.buffer() >> a;
      alpha[np] = std::min( alpha[np], a ); //find smallest alpha between neighbors
   }

   // calculate flow for every edge (process-process connection) for every level
   
   std::map< mpi::MPIRank, mpi::MPISize > ranksToRecvFrom;
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      ranksToRecvFrom[ static_cast< mpi::MPIRank >(*n) ] = numeric_cast< mpi::MPISize >( levelsToProcess * mpi::BufferSizeTrait<double>::size );

   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1709 ); // dynamicdiffusion = 100 121 110 097 109 105 099 100 105 102 102 117 115 105 111 110 + 1
   bufferSystem.setReceiverInfo( ranksToRecvFrom );

   std::map< uint_t, std::vector< double > > flow; //process rank -> flow on every level
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      flow[*n].resize( levels, 0.0 );
   
   std::vector< double > localWeight( processWeight ); //per level
   
   double flowIterations( double_c( flowIterations_ ) );
   if( iteration >= flowIterationsIncreaseStart_ )
      flowIterations += double_c( iteration + uint_t(1) - flowIterationsIncreaseStart_ ) * flowIterationsIncrease_;
   
   const uint_t iterations = uint_c( flowIterations + 0.5 );
   for( uint_t i = uint_t(0); i < iterations; ++i )
   {
      WALBERLA_ASSERT_EQUAL( localWeight.size(), levels );

      for( auto rank = ranksToRecvFrom.begin(); rank != ranksToRecvFrom.end(); ++rank )
         for( uint_t l = uint_t(0); l < levels; ++l )
            if( processLevel[l] )
               bufferSystem.sendBuffer( rank->first ) << localWeight[l];
      
      bufferSystem.sendAll();
      
      std::vector< double > previousLocalWeight( localWeight );
      
      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         const uint_t np = uint_c( recvIt.rank() );
         for( uint_t l = uint_t(0); l < levels; ++l )
         {
            if( processLevel[l] )
            {
               WALBERLA_ASSERT( flow.find( np ) != flow.end() );
               WALBERLA_ASSERT( alpha.find( np ) != alpha.end() );
               WALBERLA_ASSERT_LESS( l, flow[np].size() );
               WALBERLA_ASSERT_LESS( l, localWeight.size() );

               double nWeight( 0.0 );
               recvIt.buffer() >> nWeight;

               const double f = alpha[ np ] * ( previousLocalWeight[l] - nWeight );
               flow[ np ][l] += f;
               localWeight[l] -= f;
            }
         }
      }
   }
   
   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
      it->second = blockforest.getProcess();
   
   if( mode_ == DIFFUSION_PUSH || ( mode_ == DIFFUSION_PUSHPULL && (iteration & uint_t(1)) == uint_t(0) ) )  // sending processes decide which blocks are exchanged
   {      
      // calculate accumulated outflow for every level
      
      std::vector< double > outflow( levels, 0.0 );
      for( auto it = flow.begin(); it != flow.end(); ++it )
         for( uint_t l = uint_t(0); l < levels; ++l )
            outflow[l] += std::max( it->second[l], 0.0 );
      
      // use global information to adapt outflow
      
      if( adaptOutflowWithGlobalInformation_ )
      {
         std::vector< double > inflow( levels, 0.0 );
         for( auto it = flow.begin(); it != flow.end(); ++it )
            for( uint_t l = uint_t(0); l < levels; ++l )
               inflow[l] += std::min( it->second[l], 0.0 );

         std::vector< double > flowScaleFactor( levels, 1.0 );
         for( uint_t l = uint_t(0); l < levels; ++l )
         {
            //if( processLevel[l] && avgProcessWeight[l] < (0.99 * processWeightLimit[l]) )
            //{
            //   const double correctedOutflow = std::max( outflow[l] - processWeightLimit[l] + avgProcessWeight[l], 0.0 );
            //   flowScaleFactor[l] = correctedOutflow / outflow[l];
            //   outflow[l] = correctedOutflow;
            //}
            
            //const double diffusionAvgWeight = processWeight[l] - outflow[l] - inflow[l];
            if( processLevel[l] )
            {
               const double correctedOutflow = std::max( processWeight[l] - processWeightLimit[l] - inflow[l], 0.0 ); // identical to below ...
               //const double correctedOutflow = std::max( outflow[l] - processWeightLimit[l] + diffusionAvgWeight, 0.0 );
               flowScaleFactor[l] = correctedOutflow / outflow[l];
               if (std::isnan(flowScaleFactor[l]))
               {
                  flowScaleFactor[l] = real_t(1);
                  continue;
               }
               outflow[l] = correctedOutflow;
            }
         }

         for( auto it = flow.begin(); it != flow.end(); ++it )
            for( uint_t l = uint_t(0); l < levels; ++l )
               if( it->second[l] > 0.0 )
                  it->second[l] *= flowScaleFactor[l];
      }
      
      // determine which blocks are send to which process

      for( uint_t l = uint_t(0); l < levels; ++l )
      {
         outflow[l] = std::min( outflow[l], processWeight[l] );
         const double outflowExcess = ( outflowExceedFactor_ - 1.0 ) * outflow[l];

         while( processLevel[l] && outflow[l] > 0.0 )
         {
            auto it = flow.begin();
            WALBERLA_ASSERT( it != flow.end() );
            WALBERLA_ASSERT_LESS( l, it->second.size() );

            // pick process with highest flow
            
            uint_t pickedProcess( it->first );
            double edgeFlow( it->second[l] );
            for( ; it != flow.end(); ++it )
            {
               WALBERLA_ASSERT_LESS( l, it->second.size() );
               if( it->second[l] > edgeFlow )
               {
                  pickedProcess = it->first;
                  edgeFlow = it->second[l];
               }
            }

            if( edgeFlow > 0.0 )
            {
               std::vector< uint_t > candidates;
               if( regardConnectivity_ && iteration < disregardConnectivityStart_ )
               {
                  // blocks of this level, not yet assigned to another process and connected to 'pickedProcess'
                  std::vector< std::set< uint_t > > connectionType( 8 ); // 0 = full face (same level), 1 = full face (different level), 2 = part face,
                                                                         // 3 = full edge (same level), 4 = full edge (different level), 5 = part edge,
                                                                         // 6 = corner (same level), 7 = corner (different level)
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const auto * block = targetProcess[i].first;
                     const uint_t level = levelwise_ ? block->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() )
                     {
                        auto faces = blockforest::getFaceNeighborhoodSectionIndices();
                        for( uint_t j = uint_t(0); j != faces.size(); ++j )
                        {
                           const uint_t sectionIndex = faces[j];
                           const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                           for( uint_t s = uint_t(0); s != sectionSize; ++s )
                           {
                              if( block->getNeighborProcess( sectionIndex, s ) == pickedProcess )
                              {
                                 if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                                    connectionType[0].insert(i);
                                 else if( block->neighborhoodSectionHasLargerBlock( sectionIndex ) )
                                    connectionType[1].insert(i);
                                 else
                                    connectionType[2].insert(i);
                              }
                           }
                        }
                        auto edges = blockforest::getEdgeNeighborhoodSectionIndices();
                        for( uint_t j = uint_t(0); j != edges.size(); ++j )
                        {
                           const uint_t sectionIndex = edges[j];
                           const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                           for( uint_t s = uint_t(0); s != sectionSize; ++s )
                           {
                              if( block->getNeighborProcess( sectionIndex, s ) == pickedProcess )
                              {
                                 if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                                    connectionType[3].insert(i);
                                 else if( block->neighborhoodSectionHasLargerBlock( sectionIndex ) )
                                    connectionType[4].insert(i);
                                 else
                                    connectionType[5].insert(i);
                              }
                           }
                        }
                        auto corners = blockforest::getCornerNeighborhoodSectionIndices();
                        for( uint_t j = uint_t(0); j != corners.size(); ++j )
                        {
                           const uint_t sectionIndex = corners[j];
                           const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                           for( uint_t s = uint_t(0); s != sectionSize; ++s )
                           {
                              if( block->getNeighborProcess( sectionIndex, s ) == pickedProcess )
                              {
                                 if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                                    connectionType[6].insert(i);
                                 else
                                    connectionType[7].insert(i);
                              }
                           }
                        }
                     }
                  }
                  std::set< uint_t > assigned;
                  for( auto type = connectionType.begin(); type != connectionType.end(); ++type )
                  {
                     for( auto index = type->begin(); index != type->end(); ++index )
                     {
                        if( assigned.find(*index) == assigned.end() )
                        {
                           candidates.push_back(*index); // -> this order leads to a prioritization of face over edge over corner connections (if everything else is equal)
                           assigned.insert(*index);
                        }
                     }
                  }
                  /*
                  std::set< uint_t > assigned;
                  // no distinction between different connection types (face,edge,corner)
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const uint_t level = levelwise_ ? targetProcess[i].first->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() )
                     {
                        bool connectedToPickedProcess( false );
                        for( auto n = targetProcess[i].first->getNeighborhood().begin();
                                 n != targetProcess[i].first->getNeighborhood().end() && !connectedToPickedProcess; ++n )
                           connectedToPickedProcess = ( n->getProcess() == pickedProcess );
                        if( connectedToPickedProcess )
                        {
                           candidates.push_back( i );
                           assigned.insert( i );
                        }
                     }
                  }
                  */
                  // blocks of this level, not yet assigned to another process and not only connected to other blocks of this process (= no 'inner' block)
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const uint_t level = levelwise_ ? targetProcess[i].first->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() &&
                         assigned.find(i) == assigned.end() )
                     {
                        bool connectedToOtherProcesses( false );
                        for( auto n = targetProcess[i].first->getNeighborhood().begin();
                                 n != targetProcess[i].first->getNeighborhood().end() && !connectedToOtherProcesses; ++n )
                           connectedToOtherProcesses = ( n->getProcess() != blockforest.getProcess() );
                        if( connectedToOtherProcesses )
                        {
                           candidates.push_back( i );
                           assigned.insert( i );
                        }
                     }
                  }
                  // all blocks of this level, not yet assigned to another process and still assigned to this process
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const uint_t level = levelwise_ ? targetProcess[i].first->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() &&
                         assigned.find(i) == assigned.end() )
                        candidates.push_back( i );
                  }
               }
               else
               {
                  std::vector< uint_t > blocksToPick;
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const uint_t level = levelwise_ ? targetProcess[i].first->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() )
                        blocksToPick.push_back( i );
                  }
                  while( ! blocksToPick.empty() )
                  {
                     const uint_t pickedBlock = random_( uint_t(0), uint_c( blocksToPick.size() ) - uint_t(1) );
                     WALBERLA_ASSERT_LESS( pickedBlock, blocksToPick.size() );
                     candidates.push_back( blocksToPick[ pickedBlock ] );
                     blocksToPick[ pickedBlock ] = blocksToPick.back();
                     blocksToPick.pop_back();
                  }
                  /*
                  for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
                  {
                     const uint_t level = levelwise_ ? targetProcess[i].first->getLevel() : uint_t(0);
                     if( level == l && targetProcess[i].second == blockforest.getProcess() )
                        candidates.push_back( i );
                  }
                  */
               }

               // only blocks that do not exceed the entire process outflow are viable candidates
               std::vector< uint_t > viableCandidates;
               for( auto candidate = candidates.begin(); candidate != candidates.end(); ++candidate )
               {
                  if( weight( targetProcess[ *candidate ].first ) <= ( outflowExcess + outflow[l] ) ) // ( outflowExceedFactor_ * outflow[l] ) )
                  {
                     viableCandidates.push_back( *candidate );
                  }
               }

               if( ! viableCandidates.empty() )
               {
                  // the block that comes closest to the edge flow of the picked neighbor process is the final candidate for moving to this process
                  uint_t finalCandidate( uint_t(0) );
                  double diff = std::abs( edgeFlow - weight( targetProcess[ viableCandidates[0] ].first ) );
                  for( uint_t i = uint_t(1); i != viableCandidates.size(); ++i )
                  {
                     double d = std::abs( edgeFlow - weight( targetProcess[ viableCandidates[i] ].first ) );
                     if( d < diff )
                     {
                        diff = d;
                        finalCandidate = i;
                     }
                  }

                  const double w = weight( targetProcess[ viableCandidates[finalCandidate] ].first );

                  targetProcess[ viableCandidates[finalCandidate] ].second = pickedProcess;
                  flow[ pickedProcess ][l] -= w;
                  outflow[l] -= w;
               }
               else
               {
                  flow[ pickedProcess ][l] = 0.0; // nothing can be done for this edge
               }
            }
            else
            {
               outflow[l] = 0.0; // breaks the loop
            }
         }
      }
   }
   else // receiving processes decide which blocks are exchanged
   {      
      WALBERLA_ASSERT( mode_ == DIFFUSION_PULL || ( mode_ == DIFFUSION_PUSHPULL && (iteration & uint_t(1)) == uint_t(1) ) );
   
      // calculate accumulated inflow for every level
      
      std::vector< double > inflow( levels, 0.0 ); // inflow is saved as a positive number!
      for( auto it = flow.begin(); it != flow.end(); ++it )
         for( uint_t l = uint_t(0); l < levels; ++l )
            inflow[l] -= std::min( it->second[l], 0.0 );
      
      // use global information to adapt inflow
      
      if( adaptInflowWithGlobalInformation_ )
      {
         std::vector< double > outflow( levels, 0.0 );
         for( auto it = flow.begin(); it != flow.end(); ++it )
            for( uint_t l = uint_t(0); l < levels; ++l )
               outflow[l] += std::max( it->second[l], 0.0 );

         std::vector< double > flowScaleFactor( levels, 1.0 );
         for( uint_t l = uint_t(0); l < levels; ++l )
         {
            if( processLevel[l] )
            {
               const double correctedInflow = std::max( processWeightLimit[l] - processWeight[l] + outflow[l], 0.0 );           
               flowScaleFactor[l] = correctedInflow / inflow[l];
               if (std::isnan(flowScaleFactor[l]))
               {
                  flowScaleFactor[l] = real_t(1);
                  continue;
               }
               inflow[l] = correctedInflow;
            }
         }

         for( auto it = flow.begin(); it != flow.end(); ++it )
            for( uint_t l = uint_t(0); l < levels; ++l )
               if( it->second[l] < 0.0 )
                  it->second[l] *= flowScaleFactor[l];
      }
      
      // determine blocks with connections to neighbors processes, sort these blocks by their connection type, and
      // send a list of pairs <block-id,weight> to these neighbors
      
      std::map< uint_t, std::vector< std::vector< std::set< uint_t > > > > blocksForNeighbors; // sorted by level
      std::vector< std::vector< uint_t > > blocksConnectedToOtherProcesses; // sorted by level
      std::vector< std::vector< uint_t > > allBlocks; // sorted by level
      
      for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      {
         blocksForNeighbors[*n].resize( levels );
         for( uint_t l = uint_t(0); l < levels; ++l )
            blocksForNeighbors[*n][l].resize( 8 ); // 0 = full face (same level), 1 = full face (different level), 2 = part face,
                                                   // 3 = full edge (same level), 4 = full edge (different level), 5 = part edge,
                                                   // 6 = corner (same level), 7 = corner (different level)
      }
      blocksConnectedToOtherProcesses.resize( levels );
      allBlocks.resize( levels );
      
      for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
      {
         const auto * block = targetProcess[i].first;
         const uint_t level = levelwise_ ? block->getLevel() : uint_t(0);
         if( processLevel[level] )
         {
            if( regardConnectivity_ && iteration < disregardConnectivityStart_ )
            {
               auto faces = blockforest::getFaceNeighborhoodSectionIndices();
               for( uint_t j = uint_t(0); j != faces.size(); ++j )
               {
                  const uint_t sectionIndex = faces[j];
                  const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                  for( uint_t s = uint_t(0); s != sectionSize; ++s )
                  {
                     const uint_t np = block->getNeighborProcess( sectionIndex, s );
                     if( np != blockforest.getProcess() )
                     {
                        WALBERLA_ASSERT( blocksForNeighbors.find(np) != blocksForNeighbors.end() );
                        if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                           blocksForNeighbors[np][level][0].insert(i);
                        else if( block->neighborhoodSectionHasLargerBlock( sectionIndex ) )
                           blocksForNeighbors[np][level][1].insert(i);
                        else
                           blocksForNeighbors[np][level][2].insert(i);
                     }
                  }
               }
               auto edges = blockforest::getEdgeNeighborhoodSectionIndices();
               for( uint_t j = uint_t(0); j != edges.size(); ++j )
               {
                  const uint_t sectionIndex = edges[j];
                  const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                  for( uint_t s = uint_t(0); s != sectionSize; ++s )
                  {
                     const uint_t np = block->getNeighborProcess( sectionIndex, s );
                     if( np != blockforest.getProcess() )
                     {
                        WALBERLA_ASSERT( blocksForNeighbors.find(np) != blocksForNeighbors.end() );
                        if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                           blocksForNeighbors[np][level][3].insert(i);
                        else if( block->neighborhoodSectionHasLargerBlock( sectionIndex ) )
                           blocksForNeighbors[np][level][4].insert(i);
                        else
                           blocksForNeighbors[np][level][5].insert(i);
                     }
                  }
               }
               auto corners = blockforest::getCornerNeighborhoodSectionIndices();
               for( uint_t j = uint_t(0); j != corners.size(); ++j )
               {
                  const uint_t sectionIndex = corners[j];
                  const uint_t sectionSize = block->getNeighborhoodSectionSize( sectionIndex );
                  for( uint_t s = uint_t(0); s != sectionSize; ++s )
                  {
                     const uint_t np = block->getNeighborProcess( sectionIndex, s );
                     if( np != blockforest.getProcess() )
                     {
                        WALBERLA_ASSERT( blocksForNeighbors.find(np) != blocksForNeighbors.end() );
                        if( block->neighborhoodSectionHasEquallySizedBlock( sectionIndex ) )
                           blocksForNeighbors[np][level][6].insert(i);
                        else
                           blocksForNeighbors[np][level][7].insert(i);
                     }
                  }
               }
               
               bool connectedToOtherProcesses( false );
               for( auto n = block->getNeighborhood().begin(); n != block->getNeighborhood().end() && !connectedToOtherProcesses; ++n )
                  connectedToOtherProcesses = ( n->getProcess() != blockforest.getProcess() );
               if( connectedToOtherProcesses )
                  blocksConnectedToOtherProcesses[level].push_back(i);
            }
            allBlocks[level].push_back(i);
         }
      }
      
      std::map< uint_t, std::vector< std::pair< BlockID, double > > > blocksForNeighborsUnsorted;
      
      for( auto it = blocksForNeighbors.begin(); it != blocksForNeighbors.end(); ++it )
      {
         const uint_t np = it->first;
         WALBERLA_ASSERT_EQUAL( it->second.size(), levels );
         for( uint_t l = uint_t(0); l < levels; ++l )
         {
            if( processLevel[l] )
            {
               if( regardConnectivity_ && iteration < disregardConnectivityStart_ )
               {
                  std::set< uint_t > assigned;
                  for( auto type = it->second[l].begin(); type != it->second[l].end(); ++type )
                  {
                     for( auto index = type->begin(); index != type->end(); ++index )
                     {
                        if( assigned.find(*index) == assigned.end() )
                        {
                           const auto * block = targetProcess[*index].first;
                           blocksForNeighborsUnsorted[np].push_back( std::make_pair( block->getId(), weight(block) ) ); // -> this order leads to a prioritization of face over edge over corner connections
                           assigned.insert(*index);
                        }
                     }
                  }
                  for( auto index = blocksConnectedToOtherProcesses[l].begin(); index != blocksConnectedToOtherProcesses[l].end(); ++index )
                  {
                     if( assigned.find(*index) == assigned.end() )
                     {
                        const auto * block = targetProcess[*index].first;
                        blocksForNeighborsUnsorted[np].push_back( std::make_pair( block->getId(), weight(block) ) );
                        assigned.insert(*index);
                     }
                  }
                  for( auto index = allBlocks[l].begin(); index != allBlocks[l].end(); ++index )
                  {
                     if( assigned.find(*index) == assigned.end() )
                     {
                        const auto * block = targetProcess[*index].first;
                        blocksForNeighborsUnsorted[np].push_back( std::make_pair( block->getId(), weight(block) ) );
                     }
                  }
               }
               else
               {
                  std::vector< uint_t > blocksToPick;
                  for( auto index = allBlocks[l].begin(); index != allBlocks[l].end(); ++index )
                     blocksToPick.push_back( *index );
                  while( ! blocksToPick.empty() )
                  {
                     const uint_t pickedBlock = random_( uint_t(0), uint_c( blocksToPick.size() ) - uint_t(1) );
                     WALBERLA_ASSERT_LESS( pickedBlock, blocksToPick.size() );
                     const auto * block = targetProcess[ blocksToPick[ pickedBlock ] ].first;
                     blocksForNeighborsUnsorted[np].push_back( std::make_pair( block->getId(), weight(block) ) );
                     blocksToPick[ pickedBlock ] = blocksToPick.back();
                     blocksToPick.pop_back();
                  }
                  /*
                  for( auto index = allBlocks[l].begin(); index != allBlocks[l].end(); ++index )
                  {
                     const auto * block = targetProcess[*index].first;
                     blocksForNeighborsUnsorted[np].push_back( std::make_pair( block->getId(), weight(block) ) );
                  }
                  */
               }
            }
         }
      }

      std::set< mpi::MPIRank > neighborsToRecvFrom;
      for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
         neighborsToRecvFrom.insert( static_cast< mpi::MPIRank >(*n) );

      mpi::BufferSystem neighborsBufferSystem( MPIManager::instance()->comm(), 1710 ); // dynamicdiffusion = 100 121 110 097 109 105 099 100 105 102 102 117 115 105 111 110 + 2
      neighborsBufferSystem.setReceiverInfo( neighborsToRecvFrom, true );
      
      for( auto rank = neighborsToRecvFrom.begin(); rank != neighborsToRecvFrom.end(); ++rank )
         neighborsBufferSystem.sendBuffer( *rank ) << blocksForNeighborsUnsorted[ uint_c( *rank ) ];
      
      neighborsBufferSystem.sendAll();
      
      // every process receives this list from all neighbors

      std::map< uint_t, std::vector< std::pair< BlockID, double > > > blocksFromNeighborsUnsorted;
         
      for( auto recvIt = neighborsBufferSystem.begin(); recvIt != neighborsBufferSystem.end(); ++recvIt )
      {
         const uint_t np = uint_c( recvIt.rank() );
         recvIt.buffer() >> blocksFromNeighborsUnsorted[ np ];
      }
      
      std::map< uint_t, std::vector< std::vector< std::pair< BlockID, double > > > > blocksFromNeighbors; // sorted by level

      for( auto it = blocksFromNeighborsUnsorted.begin(); it != blocksFromNeighborsUnsorted.end(); ++it )
      {
         blocksFromNeighbors[ it->first ].resize( levels );
         for( auto pr = it->second.begin(); pr != it->second.end(); ++pr )
         {
            const uint_t level = levelwise_ ? blockforest.getLevelFromBlockId( pr->first ) : uint_t(0);
            WALBERLA_ASSERT_LESS( level, blocksFromNeighbors[ it->first ].size() );
            blocksFromNeighbors[ it->first ][ level ].push_back( *pr );
         }
      }
      
      // using process inflow and inflow edges -> decide which blocks we want to fetch
      
      std::map< uint_t, std::vector< BlockID > > blocksToFetch;
      std::set< BlockID > pickedBlocks;

      for( uint_t l = uint_t(0); l < levels; ++l )
      {
         const double inflowExcess = ( inflowExceedFactor_ - 1.0 ) * inflow[l];

         while( processLevel[l] && inflow[l] > 0.0 )
         {
            auto it = flow.begin();
            WALBERLA_ASSERT( it != flow.end() );
            WALBERLA_ASSERT_LESS( l, it->second.size() );

            // pick process with highest inflow
            
            uint_t pickedProcess( it->first );
            double edgeFlow( it->second[l] );
            for( ; it != flow.end(); ++it )
            {
               WALBERLA_ASSERT_LESS( l, it->second.size() );
               if( it->second[l] < edgeFlow )
               {
                  pickedProcess = it->first;
                  edgeFlow = it->second[l];
               }
            }

            if( edgeFlow < 0.0 )
            {
               WALBERLA_ASSERT( blocksFromNeighbors.find(pickedProcess) != blocksFromNeighbors.end() );
               const auto & candidates = blocksFromNeighbors[pickedProcess][l];
      
               // only blocks that do not exceed the entire process inflow are viable candidates
               std::vector< std::pair< BlockID, double > > viableCandidates;
               for( auto candidate = candidates.begin(); candidate != candidates.end(); ++candidate )
               {
                  if( pickedBlocks.find( candidate->first ) == pickedBlocks.end() && candidate->second <= ( inflowExcess + inflow[l] ) )
                  {
                     viableCandidates.push_back( *candidate );
                  }
               }

               if( ! viableCandidates.empty() )
               {
                  // the block that comes closest to the edge flow of the picked neighbor process is the final candidate for moving to this process
                  uint_t finalCandidate( uint_t(0) );
                  double diff = std::abs( edgeFlow + viableCandidates[finalCandidate].second );
                  for( uint_t i = uint_t(1); i != viableCandidates.size(); ++i )
                  {
                     double d = std::abs( edgeFlow + viableCandidates[i].second );
                     if( d < diff )
                     {
                        diff = d;
                        finalCandidate = i;
                     }
                  }

                  const double w = viableCandidates[finalCandidate].second;

                  blocksToFetch[pickedProcess].push_back( viableCandidates[finalCandidate].first );
                  pickedBlocks.insert( viableCandidates[finalCandidate].first );
                  
                  flow[ pickedProcess ][l] += w;
                  inflow[l] -= w;
               }
               else
               {
                  flow[ pickedProcess ][l] = 0.0; // nothing can be done for this edge
               }
            }
            else
            {
               inflow[l] = 0.0; // breaks the loop
            }
         }
      }
      
      // tell neighbors about blocks we want to fetch
      
      for( auto rank = neighborsToRecvFrom.begin(); rank != neighborsToRecvFrom.end(); ++rank )
         neighborsBufferSystem.sendBuffer( *rank ) << blocksToFetch[ uint_c( *rank ) ];
      
      neighborsBufferSystem.sendAll();
      
      std::map< uint_t, std::vector< BlockID > > blocksToSend;
      
      for( auto recvIt = neighborsBufferSystem.begin(); recvIt != neighborsBufferSystem.end(); ++recvIt )
      {
         const uint_t np = uint_c( recvIt.rank() );
         recvIt.buffer() >> blocksToSend[ np ];
      }
      
      // every process now knows which blocks must be sent to neighbors -> mark them (more than one neighbor might want the same block -> sort by edge outflow, highest outflow wins)
      
      std::map< BlockID, std::vector< uint_t > > processToSendTo;
      
      for( auto it = blocksToSend.begin(); it != blocksToSend.end(); ++it )
      {
         const uint_t p = it->first;
         for( auto id = it->second.begin(); id != it->second.end(); ++id )
            processToSendTo[ *id ].push_back( p );
      }
      
      for( auto it = processToSendTo.begin(); it != processToSendTo.end(); ++it )
      {
         uint_t index( uint_t(0) );
         while( index < targetProcess.size() && targetProcess[ index ].first->getId() != it->first )
         {
            ++index;
         }
         WALBERLA_ASSERT( index < targetProcess.size() && targetProcess[ index ].first->getId() == it->first );
         const uint_t level = levelwise_ ? targetProcess[ index ].first->getLevel() : uint_t(0);
         
         auto p = it->second.begin();
         uint_t pickedProcess = *p;
         for( ++p; p != it->second.end(); ++p )
         {
            WALBERLA_ASSERT( flow.find(*p) != flow.end() );
            if( flow[*p][level] > flow[pickedProcess][level] )
               pickedProcess = *p;
         }
         
         targetProcess[ index ].second = pickedProcess;
      }
   }

   // synchronize with neighbors if blocks are exchanged

   std::map< uint_t, uint8_t > sendBlocksToNeighbor;
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      sendBlocksToNeighbor[*n] = uint8_t(0);

   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
      if( it->second != blockforest.getProcess() )
         sendBlocksToNeighbor[ it->second ] = uint8_t(1);

   std::map< mpi::MPIRank, mpi::MPISize > blocksRanksToRecvFrom;
   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      blocksRanksToRecvFrom[ static_cast< mpi::MPIRank >(*n) ] = mpi::BufferSizeTrait<uint8_t>::size;

   mpi::BufferSystem blocksBufferSystem( MPIManager::instance()->comm(), 1711 ); // dynamicdiffusion = 100 121 110 097 109 105 099 100 105 102 102 117 115 105 111 110 + 3
   blocksBufferSystem.setReceiverInfo( blocksRanksToRecvFrom );

   for( auto rank = blocksRanksToRecvFrom.begin(); rank != blocksRanksToRecvFrom.end(); ++rank )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c(rank->first), blockforest.getProcess() );
      WALBERLA_ASSERT( sendBlocksToNeighbor.find( uint_c(rank->first) ) != sendBlocksToNeighbor.end() );
      blocksBufferSystem.sendBuffer( rank->first ) << sendBlocksToNeighbor[ uint_c(rank->first) ];
   }

   blocksBufferSystem.sendAll();

   for( auto recvIt = blocksBufferSystem.begin(); recvIt != blocksBufferSystem.end(); ++recvIt )
   {
      uint_t np = uint_c( recvIt.rank() );
      WALBERLA_ASSERT_UNEQUAL( np, blockforest.getProcess() );
      uint8_t boolean;
      recvIt.buffer() >> boolean;
      if( boolean == uint8_t(1) )
         processesToRecvFrom.insert(np);
   }

   return ( iteration + uint_t(1) ) < maxIterations_;
}

///This class is deprecated use DynamicDiffusionBalance instead.
template< typename PhantomData_T >
using DynamicLevelwiseDiffusionBalance [[deprecated]] = DynamicDiffusionBalance<PhantomData_T> ;

} // namespace blockforest
} // namespace walberla
