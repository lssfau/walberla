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
//! \file BoundaryFromVoxelFile.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Initializer.h"
#include "BoundarySetter.h"

#include "geometry/structured/VoxelFileReader.h"

#include "boundary/Boundary.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include "domain_decomposition/BlockStorage.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/IBlockID.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

#include <map>
#include <utility>
#include <vector>


namespace walberla {
namespace geometry {
namespace initializer {


struct IBlockIDPtrCompare {
   bool operator()( const IBlockID * lhs, const IBlockID * rhs ) const;
};

typedef std::map<const IBlockID*, CellInterval, IBlockIDPtrCompare> CellIntervalMap;
typedef std::map<const IBlockID*, std::pair<CellInterval, std::vector<uint8_t> >, IBlockIDPtrCompare> CellIntervalDataMap;


//*******************************************************************************************************************
/*! Sets boundary conditions using information obtained from a voxel file.

\verbatim
{
    file    pathToVoxelFile.dat;
    offset  <5,5,5>;

    Flag
    {
       value 2;
       <NameOfBoundary> {}

    }
}
\endverbatim
*
* \ingroup geometry
*/
//*******************************************************************************************************************
template <typename BoundaryHandlerT>
class BoundaryFromVoxelFile : public Initializer
{
public:
   BoundaryFromVoxelFile( const StructuredBlockStorage & structuredBlockStorage, BlockDataID & boundaryHandlerID );

   virtual void init( BlockStorage & blockStorage, const Config::BlockHandle & blockHandle );


protected:

   CellIntervalMap getIntersectedCellIntervals( const std::string & geometryFile, const Cell & offset ) const;

   BlockDataID                    boundaryHandlerID_;
   const StructuredBlockStorage & structuredBlockStorage_;
};




CellIntervalDataMap readCellIntervalsOnRoot( const std::string & geometryFile, const Cell & offset,
                                             const CellIntervalMap & cellIntervals );

CellVector findCellsWithFlag( const CellInterval & cellInterval, const std::vector<uint8_t> & data, uint8_t flag );


} // namespace initializer
} // namespace geometry
} // namespace walberla

#include "BoundaryFromVoxelFile.impl.h"
