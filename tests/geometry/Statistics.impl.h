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
//! \file Statistics.impl.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "Statistics.h"

#include "bc/HandlingCollectionInterface.h"
#include "bc/UIDFunctions.h"

#include "core/debug/Debug.h"
#include "core/mpi/MPIWrapper.h"

#include "lbm_bc/UIDFunctions.h"

#include <fstream>
#include <limits>
#include <ostream>


namespace walberla {

template< typename FlagFieldT >
void BlockStatistics<FlagFieldT>::clear()
{
   numCells_.clear();
   numLiquidCells_.clear();
   liquidCellFraction_.clear();
   overheadFraction_.clear();
}

template< typename FlagFieldT >
void BlockStatistics<FlagFieldT>::updateOnProcess()
{
   clear();

   for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
   {
      // Every process counts its allocated blocks stuff

      real_t cellCount = real_c( structuredBlockStorage_->getNumberOfXCells(*blockIt) * structuredBlockStorage_->getNumberOfYCells(*blockIt) * structuredBlockStorage_->getNumberOfZCells(*blockIt) );
      numCells_.insert( cellCount );

      FlagFieldT * flagField = blockIt->getData<FlagFieldT>( flagFieldID_ );
      typename FlagFieldT::flag_t liquidFlag = flagField->getFlag( liquidFlagUID_ );

      uint_t ctr = 0;
      for( auto it = flagField->begin(); it != flagField->end(); ++it )
         if( *it & liquidFlag )
            ++ctr;

      real_t liquidCellCount = real_c(ctr);
      numLiquidCells_.insert(liquidCellCount);

      liquidCellFraction_.insert(liquidCellCount / cellCount);
      overheadFraction_.insert(1-liquidCellCount / cellCount);

   }
}

template< typename FlagFieldT >
void BlockStatistics<FlagFieldT>::mpiAllGather()
{
   numCells_.mpiAllGather();
   numLiquidCells_.mpiAllGather();
   liquidCellFraction_.mpiAllGather();
   overheadFraction_.mpiAllGather();

   WALBERLA_ASSERT_EQUAL( numCells_.size(),           numLiquidCells_.size()     );
   WALBERLA_ASSERT_EQUAL( numLiquidCells_.size(),     liquidCellFraction_.size() );
   WALBERLA_ASSERT_EQUAL( liquidCellFraction_.size(), overheadFraction_.size()   );

}

template< typename FlagFieldT >
void ProcessStatistics<FlagFieldT>::clear()
{
   numAllocatedBlocks_.clear();
   numCells_.clear();
   numLiquidCells_.clear();
   liquidCellFraction_.clear();
   overheadFraction_.clear();
}

template< typename FlagFieldT >
void ProcessStatistics<FlagFieldT>::updateOnProcess()
{
   clear();

   blockStatistics_.updateOnProcess();

   numAllocatedBlocks_.insert( real_c( structuredBlockStorage_->getNumberOfBlocks() ) );

   real_t cellCount           = blockStatistics_.numCells().sum();
   real_t liquidCellCount     = blockStatistics_.numLiquidCells().sum();
   real_t fractionLiquidCells = math::equal( cellCount, real_t(0) ) ? real_c(1) : liquidCellCount / cellCount;
   WALBERLA_ASSERT_GREATER_EQUAL( cellCount, liquidCellCount );

   numCells_.insert( cellCount );
   numLiquidCells_.insert( liquidCellCount );
   liquidCellFraction_.insert( fractionLiquidCells );
   overheadFraction_.insert( real_c(1)-fractionLiquidCells );
}

template< typename FlagFieldT >
void ProcessStatistics<FlagFieldT>::mpiAllGather()
{
   numAllocatedBlocks_.mpiAllGather();
   numCells_.mpiAllGather();
   numLiquidCells_.mpiAllGather();
   liquidCellFraction_.mpiAllGather();
   overheadFraction_.mpiAllGather();

   blockStatistics_.mpiAllGather();

   WALBERLA_ASSERT_EQUAL( numAllocatedBlocks_.size(),  numCells_.size()           );
   WALBERLA_ASSERT_EQUAL( numCells_.size(),            numLiquidCells_.size()     );
   WALBERLA_ASSERT_EQUAL( numLiquidCells_.size(),      liquidCellFraction_.size() );
   WALBERLA_ASSERT_EQUAL( liquidCellFraction_.size(),  overheadFraction_.size()   );
}

template< typename FlagFieldT >
DomainStatistics<FlagFieldT>::DomainStatistics( const shared_ptr<StructuredBlockForest> & structuredBlockForest, BlockDataID flagFieldID, FlagUID liquidFlagUID )
   : numProcesses_(0), numBlocks_(0), numAllocatedBlocks_(0), allocatedBlocksFraction_(0), numCells_(0),
     numLiquidCells_(0), numCellsAllocated_(0), liquidCellFraction_(0), allocatedCellFraction_(0), liquidCellFromAllocatedFraction_(0), overheadFraction_(0),
     processStatistics_( structuredBlockForest, flagFieldID, liquidFlagUID ), structuredBlockForest_( structuredBlockForest )
{

}

template< typename FlagFieldT >
void DomainStatistics<FlagFieldT>::update()
{
   processStatistics_.updateOnProcess();
   processStatistics_.mpiAllGather();

   numProcesses_ = real_c( MPIManager::instance()->numProcesses() );

   typedef Vector3<real_t> Vec3;

   blockFieldSize_     = Vec3( real_c( structuredBlockForest_->getXSize() ), real_c( structuredBlockForest_->getYSize() ), real_c( structuredBlockForest_->getZSize() ) );
   numBlocks_          = blockFieldSize_[0] * blockFieldSize_[1] * blockFieldSize_[2];
   blockSize_          = Vec3( real_c( structuredBlockForest_->getNumberOfXCells( *structuredBlockForest_->begin() ) ),
                               real_c( structuredBlockForest_->getNumberOfYCells( *structuredBlockForest_->begin() ) ),
                               real_c( structuredBlockForest_->getNumberOfZCells( *structuredBlockForest_->begin() ) ) );
   size_               = Vec3( blockFieldSize_[0] * blockSize_[0], blockFieldSize_[1] * blockSize_[1], blockFieldSize_[2] * blockSize_[2] );

   numAllocatedBlocks_ = processStatistics_.numAllocatedBlocks().sum();
   WALBERLA_ASSERT_GREATER      ( blockFieldSize_[0], 0 );
   WALBERLA_ASSERT_GREATER      ( blockFieldSize_[1], 0 );
   WALBERLA_ASSERT_GREATER      ( blockFieldSize_[2], 0 );
   WALBERLA_ASSERT_GREATER      ( numBlocks_, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( numBlocks_, numAllocatedBlocks_ );
   allocatedBlocksFraction_ = numAllocatedBlocks_ / numBlocks_;

   numCells_          = size_[0] * size_[1] * size_[2];
   numLiquidCells_    = processStatistics_.numLiquidCells().sum();
   numCellsAllocated_ = processStatistics_.numCells().sum();

   WALBERLA_ASSERT_GREATER      ( numCells_,          0                  );
   WALBERLA_ASSERT_GREATER_EQUAL( numCells_,          numCellsAllocated_ );
   WALBERLA_ASSERT_GREATER_EQUAL( numCellsAllocated_, numLiquidCells_    );

   liquidCellFraction_              = numLiquidCells_ / numCells_;
   allocatedCellFraction_           = numCellsAllocated_ / numCells_;
   liquidCellFromAllocatedFraction_ = numLiquidCells_ / numCellsAllocated_;
   overheadFraction_                = real_c(1) - liquidCellFromAllocatedFraction_;
}

template< typename FlagFieldT >
std::ostream & operator<<( std::ostream & os, const DomainStatistics<FlagFieldT> & ds)
{
   os << "Domain Statistics:\n"
      << "==================\n"
      << "Blockfield size:              " << ds.blockFieldSize_          << "\n"
      << "Block size:                   " << ds.blockSize_               << "\n"
      << "Number of blocks:             " << ds.numBlocks_               << "\n"
      << "Number of allocated blocks:   " << ds.numAllocatedBlocks_      << "\n"
      << "Fraction of allocated blocks: " << ds.allocatedBlocksFraction_ << "\n"
      << "\n"
      << "Domain size:               " << ds.size_              << "\n"
      << "Number of cells:           " << ds.numCells_          << "\n"
      << "Number of liquid cells:    " << ds.numLiquidCells_    << "\n"
      << "Number of allocated cells: " << ds.numCellsAllocated_ << "\n"
      << "\n"
      << "Fraction of liquid cells:    " << ds.liquidCellFraction_    << "\n"
      << "Fraction of allocated cells: " << ds.allocatedCellFraction_ << "\n"
      << "\n"
      << "Fraction of liquid cells from allocated cells:     " << ds.liquidCellFromAllocatedFraction_ << "\n"
      << "Fraction of non-liquid cells from allocated cells: " << ds.overheadFraction_                << "\n"
      << "\n";

   os << ds.processStatistics_;

   return os;
}

template< typename FlagFieldT >
std::ostream & operator<<(std::ostream & os, const ProcessStatistics<FlagFieldT> & ps)
{
   static const std::string formatString("avg: %mean stdDev: %stddev min: %min max: %max");

   os << "Process Statistics:\n"
      << "==================\n"
      << "Number of allocated blocks: " << ps.numAllocatedBlocks().format(formatString) << "\n"
      << "\n"
      << "Number of cells:        " << ps.numCells().format(formatString)       << "\n"
      << "Number of liquid cells: " << ps.numLiquidCells().format(formatString) << "\n"
      << "\n"
      << "Fraction of liquid cells:   " << ps.liquidCellFraction().format(formatString) << "\n"
      << "Fraction of overhead cells: " << ps.overheadFraction().format(formatString)   << "\n"
      << "\n";

   os << ps.blockStatistics();

   return os;
}

template< typename FlagFieldT >
std::ostream & operator<<(std::ostream & os, const BlockStatistics<FlagFieldT> & bs)
{
   static const std::string formatString("avg: %mean stdDev: %stddev min: %min max: %max");

   os << "Allocated Block Statistics:\n"
      << "===========================\n"
      << "Number of cells:           " << bs.numCells().format(formatString)          << "\n"
      << "Number of liquid cells:    " << bs.numLiquidCells().format(formatString)    << "\n"
      << "\n"
      << "Fraction of liquid cells:    " << bs.liquidCellFraction().format(formatString)    << "\n"
      << "Fraction of overhead cells:  " << bs.overheadFraction().format(formatString)      << "\n";

   return os;
}

void writeToStreamAsLine( std::ostream & os, const math::Sample & statReal, const std::string & name)
{
   os << "# " << name << "\n";
   std::copy( statReal.begin(), statReal.end(), std::ostream_iterator<real_t>(os, " ") );
   os << "\n";
}

template< typename FlagFieldT >
void writeToFile(const std::string & filename, const ProcessStatistics<FlagFieldT> & ps)
{
   std::ofstream os(filename.c_str());
   writeToStreamAsLine(os, ps.numAllocatedBlocks(), "numAllocatedBlocks");
   writeToStreamAsLine(os, ps.numCells(), "numCells");
   writeToStreamAsLine(os, ps.numLiquidCells(), "numLiquidCells");
   writeToStreamAsLine(os, ps.liquidCellFraction(), "liquidCellFraction");
   writeToStreamAsLine(os, ps.overheadFraction(), "overheadFraction");
}

template< typename FlagFieldT >
void writeToFile(const std::string & filename, const BlockStatistics<FlagFieldT> & bs)
{
   std::ofstream os(filename.c_str());
   writeToStreamAsLine(os, bs.numCells(), "numCells");
   writeToStreamAsLine(os, bs.numLiquidCells(), "numLiquidCells");
   writeToStreamAsLine(os, bs.liquidCellFraction(), "liquidCellFraction");
   writeToStreamAsLine(os, bs.overheadFraction(), "overheadFraction");
}

template< typename FlagFieldT >
void writeToFile(const std::string & filename, const DomainStatistics<FlagFieldT> & ds)
{
   std::ofstream os(filename.c_str());
   os << "blockFieldSizeX " << ds.blockFieldSize_[0] << "\n"
      << "blockFieldSizeY " << ds.blockFieldSize_[1] << "\n"
      << "blockFieldSizeZ " << ds.blockFieldSize_[2] << "\n"
      << "numBlocks       " << ds.numBlocks_         << "\n"
      << "numProccesses   " << ds.numProcesses_      << "\n";
}

} // namespace walberla
