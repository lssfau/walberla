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
//! \file Statistics.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/math/Sample.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

#include <iosfwd>


namespace walberla {

template< typename FlagFieldT >
class BlockStatistics
{
public:
	BlockStatistics( const shared_ptr<StructuredBlockStorage> & structuredBlockStorage, BlockDataID flagFieldID, FlagUID liquidFlagUID )
		: structuredBlockStorage_( structuredBlockStorage ), flagFieldID_(flagFieldID), liquidFlagUID_(liquidFlagUID) { }

	void updateOnProcess();
	void mpiAllGather();
	void clear();

	const math::Sample & numCells()           const { return numCells_; }
	const math::Sample & numLiquidCells()     const { return numLiquidCells_; }
	const math::Sample & liquidCellFraction() const { return liquidCellFraction_; }
	const math::Sample & overheadFraction()   const { return overheadFraction_; }

protected:
	math::Sample  numCells_;
	math::Sample  numLiquidCells_;
	math::Sample  liquidCellFraction_;
	math::Sample  overheadFraction_;

	shared_ptr<StructuredBlockStorage> structuredBlockStorage_;
	BlockDataID flagFieldID_;
	FlagUID     liquidFlagUID_;
};

template< typename FlagFieldT >
class ProcessStatistics
{
public:
	ProcessStatistics( const shared_ptr<StructuredBlockStorage> & structuredBlockStorage, BlockDataID flagFieldID, FlagUID liquidFlagUID )
		: structuredBlockStorage_( structuredBlockStorage ), blockStatistics_( structuredBlockStorage, flagFieldID, liquidFlagUID ) { }

	void updateOnProcess();
	void mpiAllGather();
	void clear();

	const math::Sample & numAllocatedBlocks()      const { return numAllocatedBlocks_; }
	const math::Sample & numCells()                const { return numCells_; }
	const math::Sample & numLiquidCells()          const { return numLiquidCells_; }
	const math::Sample & liquidCellFraction()      const { return liquidCellFraction_; }
	const math::Sample & overheadFraction()        const { return overheadFraction_; }

	const BlockStatistics<FlagFieldT> & blockStatistics()  const { return blockStatistics_; }

protected:
	math::Sample numAllocatedBlocks_;
	math::Sample numCells_;
	math::Sample numLiquidCells_;
	math::Sample liquidCellFraction_;
	math::Sample overheadFraction_;

	shared_ptr<StructuredBlockStorage> structuredBlockStorage_;

	BlockStatistics<FlagFieldT> blockStatistics_;
};

template< typename FlagFieldT >
class DomainStatistics
{
public:
	DomainStatistics( const shared_ptr<StructuredBlockForest> & structuredBlockForest, BlockDataID flagFieldID, FlagUID liquidFlagUID );
	void update();

	Vector3<real_t>   size_;
	Vector3<real_t>   blockSize_;
	Vector3<real_t>   blockFieldSize_;
	real_t numProcesses_;
	real_t numBlocks_;
	real_t numAllocatedBlocks_;
	real_t allocatedBlocksFraction_;
	real_t numCells_;
	real_t numLiquidCells_;
	real_t numCellsAllocated_;
	real_t liquidCellFraction_;
	real_t allocatedCellFraction_;
	real_t liquidCellFromAllocatedFraction_;
	real_t overheadFraction_;

	ProcessStatistics<FlagFieldT> processStatistics_;

protected:
	shared_ptr<StructuredBlockForest> structuredBlockForest_;
	BlockDataID flagFieldID_;
	FlagUID     liquidFlagUID_;
};

template< typename FlagFieldT >
std::ostream & operator<<(std::ostream & os, const DomainStatistics<FlagFieldT> & ds);

template< typename FlagFieldT >
std::ostream & operator<<(std::ostream & os, const ProcessStatistics<FlagFieldT> & ps);

template< typename FlagFieldT >
std::ostream & operator<<(std::ostream & os, const BlockStatistics<FlagFieldT> & bs);

inline void writeToStreamAsLine( std::ostream & os, const math::Sample & statReal, const std::string & name);


template< typename FlagFieldT >
void writeToFile(const std::string & filename, const ProcessStatistics<FlagFieldT> & ps);

template< typename FlagFieldT >
void writeToFile(const std::string & filename, const BlockStatistics<FlagFieldT> & bs);

inline void printStatistics();

} // namespace walberla

#include "Statistics.impl.h"
