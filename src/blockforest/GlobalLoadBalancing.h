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
//! \file GlobalLoadBalancing.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h" // for macro WALBERLA_BUILD_WITH_METIS

#include "BlockID.h"
#include "Types.h"

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/load_balancing/MetisWrapper.h"
#include "core/math/KahanSummation.h"

#include <functional>

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>


namespace walberla {
namespace blockforest {



class GlobalLoadBalancing {

public:

   template< typename BLOCK >
   class MetisConfiguration {

   public:
      using CommunicationFunction = std::function<memory_t (const BLOCK *const, const BLOCK *const)>;

      MetisConfiguration( const bool _includeMetis = false, const bool _forceMetis = false, CommunicationFunction _communicationFunction = nullptr,
                          const real_t _maxUbvec = real_c(1.5), const uint_t _iterations = uint_c(10) ) :
         includeMetis_( _includeMetis ), forceMetis_( _forceMetis ), communicationFunction_( _communicationFunction ),
         maxUbvec_( _maxUbvec ), iterations_( _iterations ) {}

      bool                  includeMetis()          const { return includeMetis_; }
      bool                  forceMetis()            const { return forceMetis_; }
      CommunicationFunction communicationFunction() const { return communicationFunction_; }
      real_t                maxUbvec()              const { return maxUbvec_; }
      uint_t                iterations()            const { return iterations_; }

   private:

      const bool                  includeMetis_;
      const bool                  forceMetis_;
            CommunicationFunction communicationFunction_;
      const real_t                maxUbvec_;
      const uint_t                iterations_;

   }; // class MetisConfiguration



   // average memory utilization (-> leads to a specific number of processes!)

   template< typename BLOCK >
   static inline uint_t balance      ( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                       const MetisConfiguration<BLOCK>& metisConfig, const memory_t avgMemoryUtilization );
   template< typename BLOCK >
   static inline uint_t balanceSorted( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                       const MetisConfiguration<BLOCK>& metisConfig, const memory_t avgMemoryUtilization );

   // number of processes

   template< typename BLOCK >
   static uint_t balance      ( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                const MetisConfiguration<BLOCK>& metisConfig, const uint_t numberOfProcesses );
   template< typename BLOCK >
   static uint_t balanceSorted( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                const MetisConfiguration<BLOCK>& metisConfig, const uint_t numberOfProcesses );

   // minimize number of processes == maximize average memory utilization

   template< typename BLOCK >
   static inline uint_t minimizeProcesses( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit,
                                           const MetisConfiguration<BLOCK>& metisConfig,
                                           const std::vector< workload_t >* processesWork = nullptr,
                                           const std::vector< memory_t >* processesMemory = nullptr );
   template< typename BLOCK >
   static inline uint_t maximizeMemoryUtilization( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit,
                                                   const MetisConfiguration<BLOCK>& metisConfig,
                                                   const std::vector< workload_t >* processesWork = nullptr,
                                                   const std::vector< memory_t >* processesMemory = nullptr );
   // optimize workload

   template< typename BLOCK >
   static uint_t optimizeWorkload( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                   const MetisConfiguration<BLOCK>& metisConfig, const real_t workloadWeighting, const bool sortByLevel = false,
                                   const uint_t samples = 10 );

   // process reordering

   template< typename BLOCK >
   static void prepareProcessReordering( const std::vector< BLOCK* > & blocks,       std::vector< std::vector< uint_t > > & processNeighbors );
   template< typename BLOCK >
   static void reorderProcessesByBFS   (       std::vector< BLOCK* > & blocks, const std::vector< std::vector< uint_t > > & processNeighbors );

private:

   // space-filling curve helper functions

   template< typename BLOCK >
   static uint_t fixedWork( const std::vector< BLOCK* >& blocks, const workload_t workloadLimit, const memory_t memoryLimit,
                            const std::vector< workload_t >* processesWork = nullptr, const std::vector< memory_t >* processesMemory = nullptr );

#ifdef WALBERLA_BUILD_WITH_METIS

   // METIS

   template< typename BLOCK >
   static uint_t metis( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit, const MetisConfiguration<BLOCK>& metisConfig,
                        const uint_t numberOfProcesses );
public:
   template< typename BLOCK, typename CommFunction >
   static void metis2( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses, const CommFunction & commFunction );
private:

   static inline uint_t metisAdaptPartVector( std::vector< int64_t >& part, const uint_t numberOfProcesses );

   template< typename BLOCK >
   static memory_t metisMaxMemory( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses, const std::vector< int64_t >& part );

   static inline std::ostream &  metisErrorCodeToStream( int errorCode, std::ostream & oss );
   static inline std::string     metisErrorCodeToString( int errorCode );

#endif // WALBERLA_BUILD_WITH_METIS

   // optimize workload helper function

   template< typename BLOCK >
   static void checkForBetterWorkloadDistributation( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses, const real_t sumWorkload,
                                                     const real_t memoryUtilization, const real_t workloadWeighting, const real_t memoryWeighting,
                                                     real_t& bestWeightedValue, std::vector< uint_t >& bestProcessMapping,
                                                     uint_t& bestNumberOfProcesses );
};



template< typename BLOCK >
inline uint_t GlobalLoadBalancing::balance( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                            const MetisConfiguration<BLOCK>& metisConfig, const memory_t avgMemoryUtilization ) {

   WALBERLA_ASSERT_LESS( static_cast< memory_t >(0), avgMemoryUtilization );
   WALBERLA_ASSERT_LESS_EQUAL( avgMemoryUtilization, static_cast< memory_t >(1) );
   WALBERLA_ASSERT_GREATER( memoryLimit, static_cast< memory_t >(0) );

   return balance( blocks, sfcIterations, memoryLimit, metisConfig,
                   uint_c( static_cast< memory_t >(0.5) + memorySum( blocks ) / ( avgMemoryUtilization * memoryLimit ) ) );
}



template< typename BLOCK >
inline uint_t GlobalLoadBalancing::balanceSorted( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                                  const MetisConfiguration<BLOCK>& metisConfig, const memory_t avgMemoryUtilization ) {

   WALBERLA_ASSERT_LESS( static_cast< memory_t >(0), avgMemoryUtilization );
   WALBERLA_ASSERT_LESS_EQUAL( avgMemoryUtilization, static_cast< memory_t >(1) );
   WALBERLA_ASSERT_GREATER( memoryLimit, static_cast< memory_t >(0) );

   return balanceSorted( blocks, sfcIterations, memoryLimit, metisConfig,
                         uint_c( static_cast< memory_t >(0.5) + memorySum( blocks ) / ( avgMemoryUtilization * memoryLimit ) ) );
}



template< typename BLOCK >
uint_t GlobalLoadBalancing::balance( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                     const MetisConfiguration<BLOCK>& metisConfig, const uint_t numberOfProcesses ) {

   WALBERLA_ASSERT_GREATER( numberOfProcesses, 0 );
   WALBERLA_ASSERT_GREATER( memoryLimit, static_cast< memory_t >(0) );

   // min workload per process

   workload_t minWorkload = workloadSum( blocks ) / static_cast< workload_t >( numberOfProcesses );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      minWorkload = std::max( minWorkload, blocks[i]->getWorkload() );
      WALBERLA_ASSERT_LESS_EQUAL( blocks[i]->getMemory(), memoryLimit );
   }

   // check min - potentially nothing more to do -> finished!

   uint_t processes = fixedWork( blocks, minWorkload, memoryLimit );

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "--------------------\nglobalLoadBalancing: static load balancing (" << blocks.size() << " blocks) to " << numberOfProcesses
   //              << " processes:\n" << "   (min) workload limit: " << minWorkload << " -> " << processes << " processes ";
   // }

#ifdef WALBERLA_BUILD_WITH_METIS
   if( processes <= numberOfProcesses && !metisConfig.forceMetis() ) {
#else
   if( processes <= numberOfProcesses ) {
#endif
      // DEBUG_LOGGING_SECTION { std::cout << "(+)\n"; }
      return processes;
   }
   // else DEBUG_LOGGING_SECTION { std::cout << "(-)\n"; }

#ifdef WALBERLA_BUILD_WITH_METIS
   if( processes > numberOfProcesses ) {
#endif

   // calculated min number of processes (respect memory limit!) -> leads to max workload

   processes = minimizeProcesses( blocks, memoryLimit, metisConfig );
   if( processes > numberOfProcesses ) {
      // DEBUG_LOGGING_SECTION { std::cout << "ERROR: static load balancing requires at least " << processes << " processes [a distribution to just "
      //                                   << numberOfProcesses << " processes is impossible]" << std::endl; }
      return 0;
   }

#ifdef WALBERLA_BUILD_WITH_METIS
   }
#endif

   workload_t maxWorkload = static_cast< workload_t >( 0 );

   std::vector< workload_t > workload( processes, static_cast< workload_t >(0) );
   for( uint_t i = 0; i != blocks.size(); ++i ) {
      WALBERLA_ASSERT_LESS( blocks[i]->getTargetProcess(), processes );
      workload[ blocks[i]->getTargetProcess() ] += blocks[i]->getWorkload();
   }
   for( uint_t i = 0; i != processes; ++i )
      maxWorkload = std::max( maxWorkload, workload[i] );

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "   (max) workload limit: " << maxWorkload << " -> " << processes << " processes (+)\n";
   // }

   // remember block<->process association

   std::vector< uint_t > blockProcessMap( blocks.size(), 0 );
   for( uint_t i = 0; i != blocks.size(); ++i )
      blockProcessMap[i] = blocks[i]->getTargetProcess();

   // binary search ...

   for( uint_t i = 0; i != sfcIterations; ++i ) {

      workload_t workloadLimit = ( minWorkload + maxWorkload ) / static_cast< workload_t >( 2 );

      uint_t p = fixedWork( blocks, workloadLimit, memoryLimit );

      // DEBUG_LOGGING_SECTION {
      //    std::cout << "         workload limit: " << workloadLimit << " -> " << p << " processes ";
      // }

      if( p <= numberOfProcesses ) {

         processes = p;

         for( uint_t j = 0; j != blocks.size(); ++j )
            blockProcessMap[j] = blocks[j]->getTargetProcess();

         maxWorkload = workloadLimit;

         // DEBUG_LOGGING_SECTION { std::cout << "(+)\n"; }
      }
      else {
         minWorkload = workloadLimit;

         // DEBUG_LOGGING_SECTION { std::cout << "(-)\n"; }
      }
   }

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      WALBERLA_ASSERT_LESS( blockProcessMap[i], numberOfProcesses );
      blocks[i]->assignTargetProcess( blockProcessMap[i] );
   }

#ifdef WALBERLA_BUILD_WITH_METIS

   if( numberOfProcesses > 1 && metisConfig.includeMetis() ) {

      maxWorkload = static_cast< workload_t >( 0 );

      workload.assign( processes, static_cast< workload_t >(0) );
      for( uint_t i = 0; i != blocks.size(); ++i ) {
         WALBERLA_ASSERT_LESS( blockProcessMap[i], processes );
         workload[ blockProcessMap[i] ] += blocks[i]->getWorkload();
      }
      for( uint_t i = 0; i != processes; ++i )
         maxWorkload = std::max( maxWorkload, workload[i] );

      uint_t p = metis( blocks, memoryLimit, metisConfig, numberOfProcesses );

      if( p != 0 ) {

         workload.assign( p, static_cast< workload_t >(0) );
         for( uint_t i = 0; i != blocks.size(); ++i ) {
            WALBERLA_ASSERT_LESS( blocks[i]->getTargetProcess(), p );
            workload[ blocks[i]->getTargetProcess() ] += blocks[i]->getWorkload();
         }
         workload_t metisMaxWorkload = static_cast< workload_t >( 0 );
         for( uint_t i = 0; i != p; ++i )
            metisMaxWorkload = std::max( metisMaxWorkload, workload[i] );

         if( metisConfig.forceMetis() || metisMaxWorkload < maxWorkload )
            processes = p;
         else {
            for( uint_t i = 0; i != blocks.size(); ++i )
               blocks[i]->assignTargetProcess( blockProcessMap[i] );
         }
      }
   }

#endif // WALBERLA_BUILD_WITH_METIS

   WALBERLA_ASSERT_LESS_EQUAL( processes, numberOfProcesses );

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "final load balancing result: blocks are distributed to " << processes << " processes [initial request: "
   //              << numberOfProcesses << ", unused processes: "
   //              << ( real_c(100) * real_c( numberOfProcesses - processes ) / real_c( numberOfProcesses ) ) << " %]\n--------------------\n";
   // }

   return processes;
}



template< typename BLOCK >
uint_t GlobalLoadBalancing::balanceSorted( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                           const MetisConfiguration<BLOCK>& metisConfig, const uint_t numberOfProcesses ) {

   WALBERLA_ASSERT_GREATER( numberOfProcesses, 0 );
   WALBERLA_ASSERT_GREATER( memoryLimit, static_cast< memory_t >(0) );

   uint_t maxLevel = uint_c(0);
   for( uint_t i = 0; i != blocks.size(); ++i )
      maxLevel = std::max( maxLevel, blocks[i]->getLevel() );

   std::vector< std::vector< BLOCK* > > sortedBlocks( maxLevel + 1 );
   for( uint_t i = 0; i != blocks.size(); ++i )
      sortedBlocks[ maxLevel - blocks[i]->getLevel() ].push_back( blocks[i] );

   std::vector< workload_t > processesWork( numberOfProcesses, static_cast< workload_t >(0) );
   std::vector< memory_t >   processesMemory( numberOfProcesses, static_cast< memory_t >(0) );

   uint_t requiredProcesses = uint_c(0);

   // WALBERLA_LOG_DETAIL_ON_ROOT( "--------------------\nglobalLoadBalancing: static load balancing (" << blocks.size() << " blocks) to "
   //                              << numberOfProcesses << " processes:" );

   for( uint_t l = 0; l != sortedBlocks.size(); ++l ) {

      // min workload per process

      workload_t minWorkload = ( workloadSum( sortedBlocks[l] ) + math::kahanSummation( processesWork.begin(), processesWork.end() ) ) /
                               static_cast< workload_t >( numberOfProcesses );

      for( uint_t i = 0; i != sortedBlocks[l].size(); ++i ) {
         minWorkload = std::max( minWorkload, sortedBlocks[l][i]->getWorkload() );
         WALBERLA_ASSERT_LESS_EQUAL( sortedBlocks[l][i]->getMemory(), memoryLimit );
      }
      for( uint_t i = 0; i != processesWork.size(); ++i )
         minWorkload = std::max( minWorkload, processesWork[i] );

      // check min - potentially nothing more to do

      uint_t processes = fixedWork( sortedBlocks[l], minWorkload, memoryLimit, &processesWork, &processesMemory );

      // WALBERLA_LOG_DETAIL_ON_ROOT( "   (min) workload limit: " << minWorkload << " -> " << processes << " processes " );

      // if( processes <= numberOfProcesses ) {
      //    DEBUG_LOGGING_SECTION { std::cout << "(+)\n"; }
      // }

      if( processes > numberOfProcesses ) {

         // calculated min number of processes (respect memory limit!) -> leads to max workload

         processes = minimizeProcesses( sortedBlocks[l], memoryLimit, metisConfig, &processesWork, &processesMemory );
         if( processes > numberOfProcesses )
            return 0;

         workload_t maxWorkload = static_cast< workload_t >( 0 );

         std::vector< workload_t > workload( processes, static_cast< workload_t >(0) );
         for( uint_t i = 0; i != sortedBlocks[l].size(); ++i ) {
            WALBERLA_ASSERT_LESS( sortedBlocks[l][i]->getTargetProcess(), processes );
            workload[ sortedBlocks[l][i]->getTargetProcess() ] += sortedBlocks[l][i]->getWorkload();
         }
         for( uint_t i = 0; i != processes; ++i )
            maxWorkload = std::max( maxWorkload, workload[i] + processesWork[i] );

         // WALBERLA_LOG_DETAIL_ON_ROOT( "   (max) workload limit: " << maxWorkload << " -> " << processes << " processes" );

         // remember block<->process association

         std::vector< uint_t > blockProcessMap( sortedBlocks[l].size(), 0 );
         for( uint_t i = 0; i != sortedBlocks[l].size(); ++i )
            blockProcessMap[i] = sortedBlocks[l][i]->getTargetProcess();

         // binary search ...

         for( uint_t i = 0; i != sfcIterations; ++i ) {

            workload_t workloadLimit = ( minWorkload + maxWorkload ) / static_cast< workload_t >( 2 );

            uint_t p = fixedWork( sortedBlocks[l], workloadLimit, memoryLimit, &processesWork, &processesMemory );

            // WALBERLA_LOG_DETAIL_ON_ROOT( "         workload limit: " << workloadLimit << " -> " << p << " processes " );

            if( p <= numberOfProcesses ) {

               processes = p;

               for( uint_t j = 0; j != sortedBlocks[l].size(); ++j )
                  blockProcessMap[j] = sortedBlocks[l][j]->getTargetProcess();

               maxWorkload = workloadLimit;
            }
            else {
               minWorkload = workloadLimit;
            }
         }

         // restore best block<->process association

         for( uint_t i = 0; i != sortedBlocks[l].size(); ++i )
            sortedBlocks[l][i]->assignTargetProcess( blockProcessMap[i] );
      }

      WALBERLA_ASSERT_LESS_EQUAL( processes, numberOfProcesses );

      for( uint_t i = 0; i != sortedBlocks[l].size(); ++i ) {
         WALBERLA_ASSERT_LESS( sortedBlocks[l][i]->getTargetProcess(), processesWork.size() );
         processesWork[ sortedBlocks[l][i]->getTargetProcess() ] += sortedBlocks[l][i]->getWorkload();
         WALBERLA_ASSERT_LESS( sortedBlocks[l][i]->getTargetProcess(), processesMemory.size() );
         processesMemory[ sortedBlocks[l][i]->getTargetProcess() ] += sortedBlocks[l][i]->getMemory();
      }

      requiredProcesses = std::max( requiredProcesses, processes );
   }

   // WALBERLA_LOG_DETAIL_ON_ROOT( "final load balancing result: blocks are distributed to " << requiredProcesses << " processes [initial request: "
   //                              << numberOfProcesses << ", unused processes: "
   //                              << ( real_c(100) * real_c( numberOfProcesses - requiredProcesses ) / real_c( numberOfProcesses ) )
   //                              << " %]\n--------------------" );

   return requiredProcesses;
}



template< typename BLOCK >
inline uint_t GlobalLoadBalancing::minimizeProcesses( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit,
#ifdef WALBERLA_BUILD_WITH_METIS
                                                      const MetisConfiguration<BLOCK>& metisConfig,
#else
                                                      const MetisConfiguration<BLOCK>& /*metisConfig*/,
#endif
                                                      const std::vector< workload_t >* processesWork,
                                                      const std::vector< memory_t >* processesMemory ) {

   // minimize number of processes == do not care about the amount of workload that is assigned to a process,
   //                                 just put as many blocks as possible on any process

   workload_t workloadLimit = workloadSum( blocks ) + ( ( processesWork == nullptr ) ? static_cast< workload_t >(0) :
                                                                                    math::kahanSummation( processesWork->begin(), processesWork->end() ) );

   uint_t numberOfProcesses = fixedWork( blocks, workloadLimit, memoryLimit, processesWork, processesMemory );

#ifdef WALBERLA_BUILD_WITH_METIS

   if( numberOfProcesses > 1 && metisConfig.includeMetis() && processesWork == NULL && processesMemory == NULL ) {

      uint_t max = numberOfProcesses;
      uint_t min = uint_c( real_c(0.5) + ( real_c( memorySum( blocks ) ) / memoryLimit ) );
      if( min == 1 ) ++min;

      if( metis( blocks, memoryLimit, metisConfig, max ) != 0 ) {
         if( metis( blocks, memoryLimit, metisConfig, min ) != 0 )
            return min;
      }

      // MAYBE-TODO: binary search in order to find optimum (i.e., minimum number of processes) ?
   }

#endif // WALBERLA_BUILD_WITH_METIS

   return numberOfProcesses;
}



template< typename BLOCK >
inline uint_t GlobalLoadBalancing::maximizeMemoryUtilization( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit,
                                                              const MetisConfiguration<BLOCK>& metisConfig,
                                                              const std::vector< workload_t >* processesWork,
                                                              const std::vector< memory_t >* processesMemory ) {

   return minimizeProcesses( blocks, memoryLimit, metisConfig, processesWork, processesMemory );
}



template< typename BLOCK >
uint_t GlobalLoadBalancing::optimizeWorkload( const std::vector< BLOCK* >& blocks, const uint_t sfcIterations, const memory_t memoryLimit,
                                              const MetisConfiguration<BLOCK>& metisConfig, const real_t workloadWeighting, const bool sortByLevel,
                                              const uint_t samples ) {

   WALBERLA_ASSERT_LESS_EQUAL( real_c(0), workloadWeighting );
   WALBERLA_ASSERT_LESS_EQUAL( workloadWeighting, real_c(1) );

   const real_t memoryWeighting = real_c(1) - workloadWeighting;

   real_t sumWorkload = real_c( workloadSum( blocks ) );
   real_t sumMemory   = real_c( memorySum( blocks ) );

   real_t                bestWeightedValue;
   std::vector< uint_t > bestProcessMapping( blocks.size(), 0 );
   uint_t                bestNumberOfProcesses = blocks.size();

   // min

   const real_t minMemoryUtilization = sumMemory / ( real_c( bestNumberOfProcesses ) * real_c( memoryLimit ) );

   real_t maxWork = real_c(0);

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      maxWork = std::max( maxWork, real_c( blocks[i]->getWorkload() ) );
      bestProcessMapping[i] = i;
   }

   const real_t avgToMaxWork = sumWorkload / ( real_c( bestNumberOfProcesses ) * maxWork );

   bestWeightedValue = avgToMaxWork * workloadWeighting + minMemoryUtilization * memoryWeighting;

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "=================== globalLoadBalancingOptimizeWorkload ===================\n=> memory utilization: "
   //              << ( real_c(100) * minMemoryUtilization ) << " %  ->  avg-work/max-work = " << avgToMaxWork
   //              << "  ->  weighted value = " << bestWeightedValue << "\n";
   //    std::cerr << ( real_c(100) * minMemoryUtilization ) << " " << ( real_c(100) * avgToMaxWork ) << " "
   //              << ( sumWorkload / real_c( bestNumberOfProcesses ) ) << " " << maxWork << "\n"; // -> gnuplot
   // }

   // max

   uint_t numberOfProcesses = maximizeMemoryUtilization( blocks, memoryLimit, metisConfig );

   const real_t maxMemoryUtilization = sumMemory / ( real_c( numberOfProcesses ) * real_c( memoryLimit ) );

   checkForBetterWorkloadDistributation( blocks, numberOfProcesses, sumWorkload, maxMemoryUtilization, workloadWeighting, memoryWeighting,
                                         bestWeightedValue, bestProcessMapping, bestNumberOfProcesses );

   // steps

   if( samples > 2 ) {

      real_t step = ( maxMemoryUtilization - minMemoryUtilization ) / real_c( samples - 1 );

      for( uint_t i = 1; i != ( samples - 1 ); ++i ) {

         const memory_t avgMemoryUtilization = static_cast< memory_t >( minMemoryUtilization + real_c(i) * step );

         if( sortByLevel )
            numberOfProcesses = balanceSorted( blocks, sfcIterations, memoryLimit, metisConfig, avgMemoryUtilization );
         else
            numberOfProcesses = balance( blocks, sfcIterations, memoryLimit, metisConfig, avgMemoryUtilization );

         checkForBetterWorkloadDistributation( blocks, numberOfProcesses, sumWorkload,
                                               sumMemory / ( real_c( numberOfProcesses ) * real_c( memoryLimit ) ), workloadWeighting, memoryWeighting,
                                               bestWeightedValue, bestProcessMapping, bestNumberOfProcesses );
      }
   }

   // restore best solution

   for( uint_t i = 0; i != blocks.size(); ++i )
      blocks[i]->assignTargetProcess( bestProcessMapping[i] );

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "\n=> chosen distribution: " << bestNumberOfProcesses << " processes (memory utilization: "
   //              << ( sumMemory / ( real_c( bestNumberOfProcesses ) * real_c( memoryLimit ) ) ) << " %, weighted value = "
   //              << bestWeightedValue << ")\n===========================================================================\n";
   // }

   return bestNumberOfProcesses;
}



template< typename BLOCK >
void GlobalLoadBalancing::prepareProcessReordering( const std::vector< BLOCK* > & blocks,
                                                    std::vector< std::vector< uint_t > > & processNeighbors ) {

   WALBERLA_ASSERT( !processNeighbors.empty() );

   const uint_t numberOfProcesses = processNeighbors.size();

   std::vector< std::map< BlockID, BLOCK* > > processBlocks( numberOfProcesses );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      WALBERLA_ASSERT_LESS( blocks[i]->getTargetProcess(), numberOfProcesses );
      processBlocks[ blocks[i]->getTargetProcess() ][ blocks[i]->getId() ] = blocks[i];
   }

   const int iNumberOfProcesses = int_c( numberOfProcesses );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int p = 0; p < iNumberOfProcesses; ++p ) {

      std::set< uint_t > neighbors;

      for( typename std::map< BlockID, BLOCK* >::iterator it = processBlocks[ uint_c(p) ].begin(); it != processBlocks[ uint_c(p) ].end(); ++it ) {
         const BLOCK* const block = (*it).second;

         for( uint_t i = 0; i != block->getNeighborhoodSize(); ++i )
            if( neighbors.insert( block->getNeighborTargetProcess(i) ).second )
               processNeighbors[ uint_c(p) ].push_back( block->getNeighborTargetProcess(i) );

//       for( uint_t n = 0; n != 26; ++n )
//          for( uint_t i = 0; i != block->getNeighborhoodSectionSize(n); ++i )
//             if( neighbors.insert( block->getNeighborTargetProcess(n,i) ).second == true )
//                processNeighbors[ uint_c(p) ].push_back( block->getNeighborTargetProcess(n,i) );
      }
   }
}



template< typename BLOCK >
void GlobalLoadBalancing::reorderProcessesByBFS( std::vector< BLOCK* > & blocks, const std::vector< std::vector< uint_t > >& processNeighbors ) {

   WALBERLA_ASSERT( !processNeighbors.empty() );

   const uint_t numberOfProcesses = processNeighbors.size();

   std::vector< bool   > processed ( numberOfProcesses, false );
   std::vector< uint_t > reorderMap( numberOfProcesses, 0     );

   uint_t process = 0;
   uint_t previousStartIndex = 0;

   for(;;)
   {
      uint_t startIndex = numberOfProcesses;
      for( uint_t i = previousStartIndex; i < numberOfProcesses; ++i )
      {
         if( !processed[i] && !processNeighbors[i].empty() )
         {
            startIndex = i;
            break;
         }
      }
      if( startIndex == numberOfProcesses )
         break;

      previousStartIndex = startIndex;

      std::list< uint_t > processList;

      processList.push_back( startIndex );
      processed[startIndex] = true;

      while( !processList.empty() ) {

         uint_t p = processList.front();
         processList.pop_front();

         reorderMap[p] = process++;

         for( uint_t i = 0; i != processNeighbors[p].size(); ++i ) {
            const uint_t neighbor = processNeighbors[p][i];

            if( !processed[neighbor] ) {
               processList.push_back( neighbor );
               processed[neighbor] = true;
            }
         }
      }

   }

#ifndef NDEBUG
   for( uint_t p = 0; p != numberOfProcesses; ++p ) {
      if( !processed[p] ) {
         WALBERLA_ASSERT_EQUAL( processNeighbors[p].size(), 0 );
         reorderMap[p] = process++;
      }
   }
   WALBERLA_ASSERT_EQUAL( process, numberOfProcesses );
   for( uint_t i = 0; i != numberOfProcesses; ++i ) {
      WALBERLA_ASSERT_LESS( reorderMap[i], numberOfProcesses );
      for( uint_t j = i+1; j != numberOfProcesses; ++j )
         WALBERLA_ASSERT_UNEQUAL( reorderMap[i], reorderMap[j] );
   }
#endif

   const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < blockssize; ++i )
      blocks[ uint_c(i) ]->assignTargetProcess( reorderMap[ blocks[ uint_c(i) ]->getTargetProcess() ] );
}



template< typename BLOCK >
uint_t GlobalLoadBalancing::fixedWork( const std::vector< BLOCK* >& blocks, const workload_t workloadLimit, const memory_t memoryLimit,
                                       const std::vector< workload_t >* processesWork, const std::vector< memory_t >* processesMemory ) {

   WALBERLA_ASSERT_GREATER( workloadLimit, static_cast< workload_t >(0) );
   WALBERLA_ASSERT_GREATER( memoryLimit  , static_cast< memory_t >(0)   );

   uint_t     processes = 0;
   workload_t workload  = ( processesWork != nullptr && processes < processesWork->size() ) ? (*processesWork)[processes] : static_cast< workload_t >(0);
   memory_t   memory    = ( processesMemory != nullptr && processes < processesMemory->size() ) ? (*processesMemory)[processes] : static_cast< memory_t >(0);

   for( uint_t i = 0; i != blocks.size(); ++i ) {

      while( ( workload + blocks[i]->getWorkload() ) > workloadLimit ||
             ( memory   + blocks[i]->getMemory()   ) > memoryLimit ) {

         WALBERLA_ASSERT_GREATER( workload, static_cast< workload_t >(0) );
         WALBERLA_ASSERT_GREATER( memory  , static_cast< memory_t >(0)   );

         ++processes;

         workload = ( processesWork != nullptr && processes < processesWork->size() ) ? (*processesWork)[processes] : static_cast< workload_t >(0);
         memory   = ( processesMemory != nullptr && processes < processesMemory->size() ) ? (*processesMemory)[processes] : static_cast< memory_t >(0);
      }

      WALBERLA_ASSERT_LESS_EQUAL( blocks[i]->getWorkload() + workload, workloadLimit );
      WALBERLA_ASSERT_LESS_EQUAL( blocks[i]->getMemory()   + memory  , memoryLimit   );

      workload += blocks[i]->getWorkload();
      memory   += blocks[i]->getMemory();

      blocks[i]->assignTargetProcess( processes );
   }

   return ( processes + 1 );
}


#ifdef WALBERLA_BUILD_WITH_METIS

template< typename BLOCK >
uint_t GlobalLoadBalancing::metis( const std::vector< BLOCK* >& blocks, const memory_t memoryLimit, const MetisConfiguration<BLOCK>& metisConfig,
                                   const uint_t numberOfProcesses )
{
   uint_t nProcesses = 0;

   // translate block forest to metis graph structure

   int64_t nvtxs = numeric_cast<int64_t>( blocks.size() ); // number of vertices in the graph
   int64_t ncon  = 2;                                    // number of balancing constraints

   for( uint_t i = 0; i != blocks.size(); ++i )
      blocks[i]->setIndex(i);

   std::vector< int64_t > xadj;                         // adjacency structure ...
   std::vector< int64_t > adjncy;                       // ... of the graph
   std::vector< int64_t > vwgt( uint_c(ncon * nvtxs) ); // weights of the vertices
   std::vector< int64_t > adjwgt;                       // weights of the edges

   xadj.push_back(0);
   for( uint_t i = 0; i != blocks.size(); ++i ) {
      const BLOCK* block = blocks[i];

      int64_t next = xadj.back() + numeric_cast< int64_t >( block->getNeighborhoodSize() );
      xadj.push_back( next );

      for( uint_t j = 0; j != block->getNeighborhoodSize(); ++j ) {
         adjncy.push_back( numeric_cast<int64_t>( block->getNeighbor(j)->getIndex() ) );
         adjwgt.push_back( !metisConfig.communicationFunction() ? 1 :
                           ( numeric_cast<int64_t>( static_cast< memory_t >(0.5) + metisConfig.communicationFunction()( block, block->getNeighbor(j) ) ) ) );
      }

      vwgt[ i * 2 ]     = numeric_cast< int64_t >( static_cast< workload_t >(0.5) + blocks[i]->getWorkload() );
      vwgt[ i * 2 + 1 ] = numeric_cast< int64_t >( static_cast<  memory_t  >(0.5) + blocks[i]->getMemory() );
   }

   int64_t nparts =  numeric_cast< int64_t >( numberOfProcesses ); // number of parts to partition the graph
   int64_t objval;

   std::vector< int64_t > part( uint_c(nvtxs) );

   real_t maxUbvec = metisConfig.maxUbvec();
   real_t ubvec[]  = { real_c(1.01), maxUbvec };

   // first call to METIS: always try to balance the workload as good as possible, but allow large imbalances concerning the memory (-> ubvec[1])

   int ret = core::METIS_PartGraphRecursive( &nvtxs, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), NULL, &(adjwgt[0]), &nparts, NULL,
                                             &(ubvec[0]), NULL /*idx t *options*/, &objval, &(part[0]) );

   // if METIS was successful AND the memory limit of each process is not violated (which is highly unlikely due to a large value for ubvec[1])
   // then the algorithm is finished

   if( ret == core::METIS_OK && metisMaxMemory( blocks, numberOfProcesses, part ) <= memoryLimit ) {

      nProcesses = metisAdaptPartVector( part, numberOfProcesses );

      for( uint_t i = 0; i != blocks.size(); ++i )
      {
         WALBERLA_ASSERT_LESS( part[i], int_c( nProcesses ) );
         WALBERLA_ASSERT_GREATER_EQUAL( part[i], 0 );
         blocks[i]->assignTargetProcess( uint_c( part[i] ) );
      }

      return nProcesses;
   }

   // second call to METIS: try to balance both workload and memory as good as possible ...

   real_t minUbvec = real_t(1);
   ubvec[1] = minUbvec;

   ret = core::METIS_PartGraphRecursive( &nvtxs, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), NULL, &(adjwgt[0]), &nparts, NULL,
                                         &(ubvec[0]), NULL /*idx t *options*/, &objval, &(part[0]) );

   // ... if this doesn't work OR if the memory limit is still violated then METIS is unable to find a valid partitioning

   if( ret != core::METIS_OK ) {

      std::string error( "METIS_ERROR" );
      if( ret == core::METIS_ERROR_INPUT ) error.assign( "METIS_ERROR_INPUT" );
      else if( ret == core::METIS_ERROR_MEMORY ) error.assign( "METIS_ERROR_MEMORY" );

      // DEBUG_LOGGING_SECTION { std::cout << "ERROR: static load balancing with METIS failed (" << error << ")" << std::endl; }
      return 0;
   }
   else if( metisMaxMemory( blocks, numberOfProcesses, part ) > memoryLimit ) {

      // DEBUG_LOGGING_SECTION { std::cout << "ERROR: static load balancing with METIS failed, the memory limit (" << memoryLimit
      //                                   << ") is violated (-> " << metisMaxMemory( blocks, numberOfProcesses, part ) << ")" << std::endl; }
      return 0;
   }

   nProcesses = metisAdaptPartVector( part, numberOfProcesses );

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      WALBERLA_ASSERT_LESS( part[i], int_c( nProcesses ) );
      WALBERLA_ASSERT_GREATER_EQUAL( part[i], 0 );
      blocks[i]->assignTargetProcess( uint_c( part[i] ) );
   }

   // binary search: iteratively call METIS and try to find a value for ubvec[1] (memory imbalance) that allows the memory imbalance to be as
   //                large as possible while at the same time the resulting partitioning does not violate the per process memory limit

   for( uint_t j = 0; j != metisConfig.iterations(); ++j ) {

      ubvec[1] = ( maxUbvec + minUbvec ) / real_c(2);

      ret = core::METIS_PartGraphRecursive( &nvtxs, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), NULL, &(adjwgt[0]), &nparts, NULL,
                                            &(ubvec[0]), NULL /*idx t *options*/, &objval, &(part[0]) );

      if( ret == core::METIS_OK && metisMaxMemory( blocks, numberOfProcesses, part ) <= memoryLimit ) {

         nProcesses = metisAdaptPartVector( part, numberOfProcesses );

         for( uint_t i = 0; i != blocks.size(); ++i )
         {
            WALBERLA_ASSERT_LESS( part[i], int_c( nProcesses ) );
            WALBERLA_ASSERT_GREATER_EQUAL( part[i], 0 );
            blocks[i]->assignTargetProcess( uint_c( part[i] ) );
         }

         minUbvec = ubvec[1];
      }
      else {
         maxUbvec = ubvec[1];
      }
   }

   return nProcesses;
}



template< typename BLOCK, typename CommFunction >
void GlobalLoadBalancing::metis2( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses, const CommFunction & commFunction )
{
   if( blocks.empty() )
      return;
  
   int64_t nvtxs = 0; // number of vertices in the graph
   int64_t ncon  = 1; // number of balancing constraints
   
   uint_t j = 0;
   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      BLOCK* block = blocks[ i ];
      if( realIsIdentical( block->getWorkload(), workload_t(0) ) )
         continue;

      block->setIndex( j++ );
   }

   std::vector< std::pair< const BLOCK*, const BLOCK* > > blockPairs;
   std::vector< real_t >                      communicationWeights;

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      const BLOCK* block = blocks[ i ];
      if( realIsIdentical( block->getWorkload(), workload_t(0) ) )
         continue;

      for( uint_t k = 0; k != blocks[ i ]->getNeighborhoodSize(); ++k )
      {
         if( block->getNeighbor( k )->getLevel() == block->getLevel() && block->getNeighbor( k )->getWorkload() > workload_t(0) )
         {
            blockPairs.push_back( std::make_pair( blocks[ i ], blocks[ i ]->getNeighbor( k ) ) );
            communicationWeights.push_back( real_t( 1 ) );
         }
      }
   }

   commFunction( blockPairs, communicationWeights );
   
   std::vector< int64_t > xadj;   // adjacency structure ...
   std::vector< int64_t > adjncy; // ... of the graph
   std::vector< int64_t > vwgt;   // weights of the vertices
   std::vector< int64_t > adjwgt; // weights of the edges
   
   xadj.push_back( 0 );
   uint_t commIdx = 0;
   for( uint_t i = 0; i != blocks.size(); ++i ) 
   {
      const BLOCK* block = blocks[ i ];

      if( realIsIdentical( block->getWorkload(), workload_t(0) ) )
         continue;
   
      ++nvtxs;
      xadj.push_back( xadj.back() );
      vwgt.push_back( numeric_cast<int64_t>( blocks[ i ]->getWorkload() ) );
   
      for( uint_t k = 0; k != block->getNeighborhoodSize(); ++k ) 
      {
         if( block->getNeighbor( k )->getLevel() == block->getLevel() && block->getNeighbor( k )->getWorkload() > workload_t(0) )
         {
            if( communicationWeights[ commIdx ] > real_t(0) )
            {
               adjncy.push_back( numeric_cast<int64_t>( block->getNeighbor( k )->getIndex() ) );
               WALBERLA_ASSERT_LESS( commIdx, communicationWeights.size() );
               WALBERLA_ASSERT_GREATER( numeric_cast<int64_t>( communicationWeights[ commIdx ] ), int64_t(0) );
               adjwgt.push_back( numeric_cast<int64_t>( communicationWeights[ commIdx ] ) );
               xadj.back() += 1;
            }
            ++commIdx;
         }
      }
   }

   WALBERLA_ASSERT_EQUAL( vwgt.size(), uint_c( nvtxs ) );
   WALBERLA_ASSERT_EQUAL( xadj.size(), uint_c( nvtxs ) + uint_t(1) );
   WALBERLA_ASSERT_EQUAL( adjncy.size(), adjwgt.size() );
   WALBERLA_ASSERT_EQUAL( adjncy.size(), communicationWeights.size() );
   
   int64_t nparts = numeric_cast<int64_t>( numberOfProcesses ); // number of parts to partition the graph
   int64_t objval;
   
   std::vector< int64_t > part( uint_c(nvtxs) );
   
   int64_t options[ METIS_NOPTIONS ];
   core::METIS_SetDefaultOptions( options );
   options[ core::METIS_OPTION_NITER ] = 1000;
   options[ core::METIS_OPTION_NSEPS ] = 100;
   options[ core::METIS_OPTION_NCUTS ] = 100;
   
   int ret = core::METIS_PartGraphKway( &nvtxs, &ncon, &( xadj[ 0 ] ), &( adjncy[ 0 ] ), &( vwgt[ 0 ] ), NULL, &( adjwgt[0] ), 
                                  &nparts, NULL, NULL, options, &objval, &( part[ 0 ] ) );

   if( ret != core::METIS_OK ) 
   {
      WALBERLA_ABORT( "METIS partitioning failed! Error: " << metisErrorCodeToString( ret ) );
   }
   
   uint_t blockCounter = 0;
   uint_t emptyBlockCounter = 0;
   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      const BLOCK* block = blocks[ i ];

      if( block->getWorkload() > workload_t(0) )
      {
         WALBERLA_ASSERT_LESS( part[ blockCounter ], int_c( numberOfProcesses ) );
         WALBERLA_ASSERT_GREATER_EQUAL( part[ blockCounter ], 0 );
         blocks[ i ]->assignTargetProcess( uint_c( part[ blockCounter++ ] ) );
      }
      else
         blocks[ i ]->assignTargetProcess( emptyBlockCounter++ % numberOfProcesses );
   }
   WALBERLA_ASSERT_EQUAL( blockCounter, uint_c( nvtxs ) );

}


std::ostream & GlobalLoadBalancing::metisErrorCodeToStream( int errorCode, std::ostream & oss )
{
   if( errorCode == core::METIS_OK)
      oss << "OK, no METIS error";
   else if( errorCode == core::METIS_ERROR_INPUT )
      oss << "Error in METIS input";
   else if (errorCode == core::METIS_ERROR_MEMORY )
      oss << "METIS could not allocate enough memory";
   else if (errorCode == core::METIS_ERROR )
      oss << "Unknown type of error";
   else
      oss << "Unknown error code";

   return oss;
}


std::string GlobalLoadBalancing::metisErrorCodeToString( int errorCode )
{
   std::ostringstream oss;
   metisErrorCodeToStream( errorCode, oss );
   return oss.str();
}



uint_t GlobalLoadBalancing::metisAdaptPartVector( std::vector< int64_t >& part, const uint_t numberOfProcesses )
{
   std::vector<bool> hasBlock( numberOfProcesses, false );
   for( uint_t i = 0; i != part.size(); ++i )
   {
      WALBERLA_ASSERT_LESS( part[i], int_c( numberOfProcesses ) );
      hasBlock[ uint_c(part[i]) ] = true;
   }

   uint_t nProcesses = 0;
   for( int i = int_c( numberOfProcesses ) - 1; i >= 0; --i )
   {
      if( !hasBlock[ uint_c(i) ] )
      {
         for( uint_t j = 0; j != part.size(); ++j )
            if( part[j] > i ) --part[j];
      }
      else
         ++nProcesses;
   }

   WALBERLA_ASSERT_LESS_EQUAL( nProcesses, numberOfProcesses );

   return nProcesses;
}



template< typename BLOCK >
memory_t GlobalLoadBalancing::metisMaxMemory( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses,
                                              const std::vector< int64_t >& part ) {

   WALBERLA_ASSERT_EQUAL( blocks.size(), part.size() );

   std::vector< memory_t > memory( numberOfProcesses, static_cast< memory_t >(0) );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      WALBERLA_ASSERT_LESS( uint_c( part[i] ), numberOfProcesses );
      memory[ uint_c( part[i] ) ] += blocks[i]->getMemory();
   }

   memory_t maxMemory = static_cast< memory_t >(0);

   for( uint_t i = 0; i != numberOfProcesses; ++i )
      maxMemory = std::max( maxMemory, memory[i] );

   return maxMemory;
}

#endif // WALBERLA_BUILD_WITH_METIS



template< typename BLOCK >
void GlobalLoadBalancing::checkForBetterWorkloadDistributation( const std::vector< BLOCK* >& blocks, const uint_t numberOfProcesses,
                                                                const real_t sumWorkload, const real_t memoryUtilization,
                                                                const real_t workloadWeighting, const real_t memoryWeighting,
                                                                real_t& bestWeightedValue, std::vector< uint_t >& bestProcessMapping,
                                                                uint_t& bestNumberOfProcesses ) {

   std::vector< real_t > workload( numberOfProcesses, real_c(0) );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      WALBERLA_ASSERT_LESS( blocks[i]->getTargetProcess(), numberOfProcesses );
      workload[ blocks[i]->getTargetProcess() ] += blocks[i]->getWorkload();
   }

   real_t maxWork = workload[0];

   for( uint_t i = 1; i != numberOfProcesses; ++i )
      maxWork = std::max( maxWork, workload[i] );

   const real_t avgToMaxWork  = sumWorkload / ( real_c( numberOfProcesses ) * maxWork );
   const real_t weightedValue = avgToMaxWork * workloadWeighting + memoryUtilization * memoryWeighting;

   // DEBUG_LOGGING_SECTION {
   //    std::cout << "=> memory utilization: " << ( real_c(100) * memoryUtilization )
   //              << " %  ->  avg-work/max-work = " << avgToMaxWork << "  ->  weighted value = " << weightedValue << "\n";
   //    std::cerr << ( real_c(100) * memoryUtilization ) << " " << ( real_c(100) * avgToMaxWork ) << " "
   //              << ( sumWorkload / real_c( numberOfProcesses ) ) << " " << maxWork << "\n"; // -> gnuplot
   // }

   if( weightedValue > bestWeightedValue ) {

      bestWeightedValue = weightedValue;

      for( uint_t i = 0; i != blocks.size(); ++i )
         bestProcessMapping[i] = blocks[i]->getTargetProcess();

      bestNumberOfProcesses = numberOfProcesses;
   }
}



} // namespace blockforest
} // namespace walberla
