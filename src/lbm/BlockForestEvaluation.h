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
//! \file BlockForestEvaluation.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/DistributedSample.h"
#include "core/mpi/MPIManager.h"
#include "core/selectable/IsSetSelected.h"
#include "core/uid/GlobalState.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/BlockCounter.h"

#include "field/CellCounter.h"
#include "field/FlagUID.h"


#include <functional>
#include <map>
#include <string>
#include <sstream>



namespace walberla {
namespace lbm {



//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the BlockForest data structure of an LBM simulation
*
*   Assumes that in-between creating an object of this class and calling any of the member functions the number of cells
*   and the number of fluid cells do not change! For simulations with static geometry, this is always the case.
*/
//**********************************************************************************************************************

template< typename CellCounter_T, typename FluidCellCounter_T, bool Pseudo2D = false >
class BlockForestEvaluationBase
{
public:

   struct BlockStatistics
   {
      real_t min;
      real_t max;
      real_t avg;
      real_t stdDev;
      real_t relStdDev;
   };

   BlockForestEvaluationBase( const weak_ptr< StructuredBlockForest > & blocks,
                              const CellCounter_T & cellCounter, const FluidCellCounter_T & fluidCellCounter,
                              const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                              const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );

   void refresh();

   void logResultOnRoot() const
   {
      WALBERLA_LOG_RESULT_ON_ROOT( "BlockForest information:\n" << loggingString() );
   }

   void logInfoOnRoot() const
   {
      WALBERLA_LOG_INFO_ON_ROOT( "BlockForest information:\n" << loggingString() );
   }

   void getResultsForSQLOnRoot( std::map< std::string, int > &         integerProperties,
                                std::map< std::string, double > &      realProperties,
                                std::map< std::string, std::string > & stringProperties );

   std::string loggingString() const;

private:
   
   //////////////////////
   // HELPER FUNCTIONS //
   //////////////////////

   static int processes() { return mpi::MPIManager::instance()->numProcesses(); }

   uint64_t allFineCells( const shared_ptr< StructuredBlockForest > & blocks ) const
   {
      uint64_t c( uint64_t(0) );
      if( Pseudo2D )
      {
         for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
            c += cells_.numberOfCells(i) * uint64_c( math::uintPow4( blocks->getNumberOfLevels() - uint_t(1) - i ) );
      }
      else
      {
         for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
            c += cells_.numberOfCells(i) * uint64_c( math::uintPow8( blocks->getNumberOfLevels() - uint_t(1) - i ) );
      }
      return c;
   }

   real_t space( const shared_ptr< StructuredBlockForest > & blocks, const uint_t level ) const
   {
      return ( Pseudo2D ? ( real_c(100) * real_c( cells_.numberOfCells(level) * uint64_c( math::uintPow4( blocks->getNumberOfLevels() - uint_t(1) - level ) ) ) ) :
                          ( real_c(100) * real_c( cells_.numberOfCells(level) * uint64_c( math::uintPow8( blocks->getNumberOfLevels() - uint_t(1) - level ) ) ) ) ) /
             real_c( allFineCells( blocks ) );
   }

   real_t memory( const uint_t level ) const
   {
      return real_c(100) * real_c( counter_.numberOfBlocks( level ) ) / real_c( counter_.numberOfBlocks() );
   }

   real_t totalWorkload( const shared_ptr< StructuredBlockForest > & blocks ) const
   {
      real_t work( real_t(0) );
      for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
         work += real_c( cells_.numberOfCells(i) * uint64_c( math::uintPow2(i) ) );
      return work;
   }

   real_t workload( const shared_ptr< StructuredBlockForest > & blocks, const uint_t level ) const
   {
      return real_c(100) * real_c( cells_.numberOfCells( level ) * uint64_c( math::uintPow2( level ) ) ) / totalWorkload( blocks );
   }

   real_t allFineWorkload( const shared_ptr< StructuredBlockForest > & blocks ) const
   {
      real_t work( real_t(0) );
      for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
         work += real_c( cells_.numberOfCells(i) * uint64_c( Pseudo2D ? math::uintPow4( blocks->getNumberOfLevels() - uint_t(1) - i ) :
                                                                        math::uintPow8( blocks->getNumberOfLevels() - uint_t(1) - i ) ) ) *
                 real_c( math::uintPow2( blocks->getNumberOfLevels() - uint_t(1) ) );
      return work;
   }



   weak_ptr< StructuredBlockForest > blocks_;

   domain_decomposition::BlockCounter counter_;
   std::vector< BlockStatistics> blockStatistics_;

   CellCounter_T cells_;
   FluidCellCounter_T fluidCells_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class BlockForestEvaluationBase

template< typename FlagField_T, bool Pseudo2D = false >
class BlockForestEvaluation : public BlockForestEvaluationBase< field::CellCounter< FlagField_T >, field::CellCounter< FlagField_T >, Pseudo2D >
{
public:
   BlockForestEvaluation( const weak_ptr< StructuredBlockForest > & blocks,
                          const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                          const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                          const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
                          : BlockForestEvaluationBase< field::CellCounter< FlagField_T >, field::CellCounter< FlagField_T >, Pseudo2D >(
                              blocks,
                              field::CellCounter< FlagField_T >( blocks, flagFieldId, Set< FlagUID >::emptySet(), requiredSelectors, incompatibleSelectors ),
                              field::CellCounter< FlagField_T >( blocks, flagFieldId, fluid, requiredSelectors, incompatibleSelectors ),
                              requiredSelectors, incompatibleSelectors )
   {
   }
};


template< typename CellCounter_T, typename FluidCellCounter_T, bool Pseudo2D >
BlockForestEvaluationBase< CellCounter_T, FluidCellCounter_T, Pseudo2D >::BlockForestEvaluationBase( const weak_ptr< StructuredBlockForest > & blocks,
                                                                                                     const CellCounter_T & cellCounter, const FluidCellCounter_T & fluidCellCounter,
                                                                                                     const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
   : blocks_( blocks ),
     counter_( blocks, requiredSelectors, incompatibleSelectors ),
     cells_( cellCounter ), fluidCells_( fluidCellCounter ),
     requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
{
   refresh();
}



template< typename CellCounter_T, typename FluidCellCounter_T, bool Pseudo2D >
void BlockForestEvaluationBase< CellCounter_T, FluidCellCounter_T, Pseudo2D >::refresh()
{
   counter_();
   cells_();
   fluidCells_();

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'BlockForestEvaluationBase' for a block storage object that doesn't exist anymore" );

   std::vector< uint_t > processBlocks( blocks->getNumberOfLevels() + uint_t(1), uint_t(0) );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const auto level = blocks->getLevel( *block );

      processBlocks[ level ] += uint_t(1);
      processBlocks.back()   += uint_t(1);
   }

   blockStatistics_.resize( blocks->getNumberOfLevels() + uint_t(1) );

   for( uint_t i = uint_t(0); i < (blocks->getNumberOfLevels() + uint_t(1)); ++i )
   {
      math::DistributedSample sample;

      sample.insert( real_c( processBlocks[i] ) );
      sample.mpiGatherRoot();

      WALBERLA_ROOT_SECTION()
      {
         blockStatistics_[i].min = sample.min();
         blockStatistics_[i].max = sample.max();
         blockStatistics_[i].avg = sample.mean();
         blockStatistics_[i].stdDev = sample.stdDeviation();
         blockStatistics_[i].relStdDev = sample.relativeStdDeviation();
      }
   }
}



template< typename CellCounter_T, typename FluidCellCounter_T, bool Pseudo2D >
void BlockForestEvaluationBase< CellCounter_T, FluidCellCounter_T, Pseudo2D >::getResultsForSQLOnRoot( std::map< std::string, int > &        integerProperties,
                                                                                                       std::map< std::string, double > &        realProperties,
                                                                                                       std::map< std::string, std::string > & stringProperties )
{
   WALBERLA_NON_ROOT_SECTION()
   {
      return;
   }

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'BlockForestEvaluationBase' for a block storage object that doesn't exist anymore" );

   realProperties[ "domainXMin" ] = blocks->getDomain().xMin();
   realProperties[ "domainXMax" ] = blocks->getDomain().xMax();
   realProperties[ "domainYMin" ] = blocks->getDomain().yMin();
   realProperties[ "domainYMax" ] = blocks->getDomain().yMax();
   realProperties[ "domainZMin" ] = blocks->getDomain().zMin();
   realProperties[ "domainZMax" ] = blocks->getDomain().zMax();

   integerProperties[ "xBlocks" ] = int_c( blocks->getXSize() );
   integerProperties[ "yBlocks" ] = int_c( blocks->getYSize() );
   integerProperties[ "zBlocks" ] = int_c( blocks->getZSize() );

   stringProperties[ "xPeriodic" ] = ( blocks->isXPeriodic() ? "yes" : "no" );
   stringProperties[ "yPeriodic" ] = ( blocks->isYPeriodic() ? "yes" : "no" );
   stringProperties[ "zPeriodic" ] = ( blocks->isZPeriodic() ? "yes" : "no" );

   integerProperties[ "levels" ] = int_c( blocks->getNumberOfLevels() );
   integerProperties[ "blocks" ] = int_c( counter_.numberOfBlocks() );

   integerProperties[ "xCells" ] = int_c( blocks->getNumberOfXCellsPerBlock() );
   integerProperties[ "yCells" ] = int_c( blocks->getNumberOfYCellsPerBlock() );
   integerProperties[ "zCells" ] = int_c( blocks->getNumberOfZCellsPerBlock() );

   realProperties[ "cells" ] = double_c( cells_.numberOfCells() );
   if( blocks->getNumberOfLevels() > uint_t(1) )
      realProperties[ "refinementCellsReduction" ] = double_c( allFineCells( blocks ) ) / real_c( cells_.numberOfCells() );
   realProperties[ "fluidCells" ] = double_c( fluidCells_.numberOfCells() );
   stringProperties[ "pseudo2D" ] = ( Pseudo2D ? "yes" : "no" );

   integerProperties[ "treeIdDigits" ]   = int_c( blocks->getTreeIdDigits() );
   integerProperties[ "blockIdBytes" ]   = int_c( blocks->getBlockIdBytes() );
   integerProperties[ "processIdBytes" ] = int_c( blocks->getProcessIdBytes() );

   realProperties[ "blocksPerProcess" ]       = blockStatistics_.back().avg;
   realProperties[ "avgBlocksPerProcess" ]    = blockStatistics_.back().avg;
   integerProperties[ "minBlocksPerProcess" ] = int_c( blockStatistics_.back().min );
   integerProperties[ "maxBlocksPerProcess" ] = int_c( blockStatistics_.back().max );

   if( blocks->getNumberOfLevels() > uint_t(1) )
   {
      for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
      {
         std::ostringstream blocksOnLevel_i;
         std::ostringstream cells_i;
         std::ostringstream fluidCells_i;
         std::ostringstream coveredSpace_i;
         std::ostringstream memoryPercentage_i;
         std::ostringstream workloadPercentage_i;
         std::ostringstream blocksPerProcess_i;
         std::ostringstream avgBlocksPerProcess_i;
         std::ostringstream minBlocksPerProcess_i;
         std::ostringstream maxBlocksPerProcess_i;

         blocksOnLevel_i << "blocksOnLevel_" << i;
         cells_i << "cells_" << i;
         fluidCells_i << "fluidCells_" << i;
         coveredSpace_i << "coveredSpace_" << i;
         memoryPercentage_i << "memoryPercentage_" << i;
         workloadPercentage_i << "workloadPercentage_" << i;
         blocksPerProcess_i << "blocksPerProcess_" << i;
         avgBlocksPerProcess_i << "avgBlocksPerProcess_" << i;
         minBlocksPerProcess_i << "minBlocksPerProcess_" << i;
         maxBlocksPerProcess_i << "maxBlocksPerProcess_" << i;

         integerProperties[ blocksOnLevel_i.str() ] = int_c( counter_.numberOfBlocks(i) );

         realProperties[ cells_i.str() ]              = double_c( cells_.numberOfCells(i) );
         realProperties[ fluidCells_i.str() ]         = double_c( fluidCells_.numberOfCells(i) );
         realProperties[ coveredSpace_i.str() ]       = double_c( space(blocks,i) );
         realProperties[ memoryPercentage_i.str() ]   = double_c( memory(i) );
         realProperties[ workloadPercentage_i.str() ] = double_c( workload(blocks,i) );

         realProperties[ blocksPerProcess_i.str() ]       = double_c( blockStatistics_[i].avg );
         realProperties[ avgBlocksPerProcess_i.str() ]    = double_c( blockStatistics_[i].avg );
         integerProperties[ minBlocksPerProcess_i.str() ] = int_c( blockStatistics_[i].min );
         integerProperties[ maxBlocksPerProcess_i.str() ] = int_c( blockStatistics_[i].max );
      }

      realProperties[ "refinementMemoryReduction" ] = real_c( allFineCells(blocks) ) / real_c( cells_.numberOfCells() );
      realProperties[ "refinementWorkloadReduction" ] = allFineWorkload(blocks) / totalWorkload(blocks);
   }
}



template< typename CellCounter_T, typename FluidCellCounter_T, bool Pseudo2D >
std::string BlockForestEvaluationBase< CellCounter_T, FluidCellCounter_T, Pseudo2D >::loggingString() const
{
   std::ostringstream oss;

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'BlockForestEvaluationBase' for a block storage object that doesn't exist anymore" );

   oss <<   "- LBM-specific block structure data:"
       << "\n   + AABB:        " << blocks->getDomain()
       << "\n   + forest size: " << blocks->getXSize() << " x " << blocks->getYSize() << " x " << blocks->getZSize() << " blocks";

   if( blocks->getNumberOfLevels() > uint_t(1) )
      oss << " (on the coarsest grid)";

   oss << "\n   + periodicity: " << ( blocks->isXPeriodic() ? "true" : "false" ) << " x "
                                 << ( blocks->isYPeriodic() ? "true" : "false" ) << " x "
                                 << ( blocks->isZPeriodic() ? "true" : "false" );

   if( blocks->getNumberOfLevels() > uint_t(1) )
      oss << "\n   + levels:      " << blocks->getNumberOfLevels();

   oss << "\n   + blocks:      " << counter_.numberOfBlocks();

   if( blocks->getNumberOfLevels() > uint_t(1) )
      oss << " (in total on all grids)";

   oss << "\n   + block size:  " << blocks->getNumberOfXCellsPerBlock() << " x " << blocks->getNumberOfYCellsPerBlock() << " x "
                                 << blocks->getNumberOfZCellsPerBlock() << " cells"
       << "\n   + cells:       " << cells_.numberOfCells();

   if( blocks->getNumberOfLevels() > uint_t(1) )
   {
      oss << " (" << allFineCells(blocks) << " if everything were fine -> data reduction by factor of "
          << ( real_c( allFineCells(blocks) ) / real_c( cells_.numberOfCells() ) ) << ")";
   }

   oss << "\n   + fluid cells: " << fluidCells_.numberOfCells() << " ("
                                 << ( real_c(100) * real_c( fluidCells_.numberOfCells() ) / real_c( cells_.numberOfCells() ) ) << " % of all cells)"
       << "\n   + pseudo 2D:   " << ( Pseudo2D ? "yes" : "no" )
       << "\n- data structure specific parameters:"
       << "\n   + tree ID digits: " << blocks->getTreeIdDigits() << " (-> block ID bytes = " << blocks->getBlockIdBytes() << ")"
       << "\n   + process ID bytes: " << blocks->getProcessIdBytes();

   if( processes() > 1 )
   {
      oss << "\n- blocks per process:"
          << "\n   + min       = " << blockStatistics_.back().min
          << "\n   + max       = " << blockStatistics_.back().max
          << "\n   + avg       = " << blockStatistics_.back().avg
          << "\n   + stdDev    = " << blockStatistics_.back().stdDev
          << "\n   + relStdDev = " << blockStatistics_.back().relStdDev;
   }

   if( blocks->getNumberOfLevels() > uint_t(1) )
   {
      oss << "\n- distribution of space/memory/work to different grid levels:";

      for( uint_t i = uint_t(0); i < blocks->getNumberOfLevels(); ++i )
      {
         oss << "\n   + level " << i
             << "\n      - " << counter_.numberOfBlocks(i) << " blocks ..."
             << "\n      - ... hold " << cells_.numberOfCells(i) << " cells"
             << "\n      - ... hold " << fluidCells_.numberOfCells(i) << " fluid cells ("
                                      << ( real_c(100) * real_c( fluidCells_.numberOfCells(i) ) / real_c( cells_.numberOfCells(i) ) )
                                      << " % of all cells on this level)"
             << "\n      - ... cover " << space(blocks,i) << " % of the total simulation space"
             << "\n      - ... account for " << memory(i) << " % of the total memory foot print"
             << "\n      - ... generate " << workload(blocks,i) << " % of the total workload *)";

         if( processes() > 1 )
         {
            oss << "\n      - blocks per process:"
                << "\n         + min       = " << blockStatistics_[i].min
                << "\n         + max       = " << blockStatistics_[i].max
                << "\n         + avg       = " << blockStatistics_[i].avg
                << "\n         + stdDev    = " << blockStatistics_[i].stdDev
                << "\n         + relStdDev = " << blockStatistics_[i].relStdDev;
         }
      }

      oss << "\n- using a uniform decomposition with a resolution equal to the finest level, one would ..."
          << "\n   + ... need " << ( real_c( allFineCells(blocks) ) / real_c( cells_.numberOfCells() ) ) << " times the memory"
          << "\n   + ... generate " << ( allFineWorkload(blocks) / totalWorkload(blocks) ) << " times the workload *)"
          << "\n\n  *) workload = total number of cell updates during the entire simulation"
          << "\n                (considers all cells, not just fluid cells)"
          << "\n                [keep in mind: on finer grids, more time steps are executed]";
   }

   return oss.str();
}



} // namespace lbm
} // namespace walberla
