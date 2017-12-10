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
//! \file BoundaryFromVoxelFile.impl.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================



namespace walberla {
namespace geometry {
namespace initializer {



template <typename BoundaryHandlerT>
BoundaryFromVoxelFile<BoundaryHandlerT>::BoundaryFromVoxelFile( const StructuredBlockStorage & structuredBlockStorage,
                                                                BlockDataID & boundaryHandlerID )
   : boundaryHandlerID_( boundaryHandlerID ), structuredBlockStorage_( structuredBlockStorage ) { }


template <typename BoundaryHandlerT>
void BoundaryFromVoxelFile<BoundaryHandlerT>::init( BlockStorage & blockStorage,
   const Config::BlockHandle & blockHandle )
{
   auto file   = blockHandle.getParameter<std::string>( "file"   );
   auto offset = blockHandle.getParameter<Cell>       ( "offset" );

   auto cellIntervals    = getIntersectedCellIntervals( file, offset );
   auto cellIntervalData = readCellIntervalsOnRoot( file, offset, cellIntervals );

   Config::Blocks configBlocks;
   blockHandle.getBlocks( configBlocks );

   std::map<uint8_t, BoundarySetter<BoundaryHandlerT> > flags;

   for( auto configBlockIt = configBlocks.begin(); configBlockIt != configBlocks.end(); ++configBlockIt )
   {
      if( configBlockIt->getKey() == "Flag" )
      {
         uint8_t flag = uint8_c( static_cast<int>( configBlockIt->template getParameter<int>( "value" ) ) );
         BoundarySetter<BoundaryHandlerT> boundarySetter;
         boundarySetter.setConfigBlock( *configBlockIt );
         flags[flag] = boundarySetter;
      }
   }

   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      auto cellIntervalDataIt = cellIntervalData.find( &(blockIt->getId()) );
      if( cellIntervalDataIt == cellIntervalData.end() )
         continue;

      for( auto setterIt = flags.begin(); setterIt != flags.end(); ++setterIt )
      {
         setterIt->second.configure( *blockIt, boundaryHandlerID_ );

         CellVector cells = findCellsWithFlag( cellIntervalDataIt->second.first,
                                               cellIntervalDataIt->second.second, setterIt->first );
         if( cells.empty() )
            continue;

         for( auto it = cells.begin(); it != cells.end(); ++it )
            structuredBlockStorage_.transformGlobalToBlockLocalCell( *it, *blockIt );

         setterIt->second.set( cells.begin(), cells.end() );
      }

   }
}


template <typename BoundaryHandlerT>
CellIntervalMap
BoundaryFromVoxelFile<BoundaryHandlerT>::getIntersectedCellIntervals( const std::string & geometryFile,
   const Cell & offset ) const
{
   // Determine global CellInterval of geometry file
   CellInterval fileCellInterval;
   WALBERLA_ROOT_SECTION()
   {
      VoxelFileReader <uint8_t> reader( geometryFile );
      fileCellInterval = CellInterval( cell_idx_c(0), cell_idx_c(0), cell_idx_c(0), cell_idx_c( reader.xSize() ) - 1,
         cell_idx_c( reader.ySize() ) - 1, cell_idx_c( reader.zSize() ) - 1 );
      fileCellInterval.shift(offset);
   }
   mpi::broadcastObject(fileCellInterval);

   // collect cell intervals of blocks intersecting the geometry file
   CellIntervalMap cellIntervals;
   for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
   {
      auto flagField = BoundarySetter<BoundaryHandlerT>::getFlagField(*blockIt, boundaryHandlerID_);
      auto ci = flagField->xyzSizeWithGhostLayer();
      structuredBlockStorage_.transformBlockLocalToGlobalCellInterval( ci, *blockIt );
      ci.intersect(fileCellInterval);
      if( !ci.empty() )
         cellIntervals[ &(blockIt->getId()) ] = ci;
   }

   return cellIntervals;
}


} // namespace initializer
} // namespace geometry
} // namespace walberla
