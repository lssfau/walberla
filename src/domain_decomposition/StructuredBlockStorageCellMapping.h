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
//! \file StructuredBlockStorageCellMapping.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StructuredBlockStorage.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace domain_decomposition {



// forward declarations

inline void transformGlobalToBlockLocal( CellVector& local, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellVector& global );
inline void transformGlobalToBlockLocal( CellVector& cells, const StructuredBlockStorage& blockStorage, const IBlock& block );

inline void transformBlockLocalToGlobal( CellVector& global, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellVector& local );
inline void transformBlockLocalToGlobal( CellVector& cells,  const StructuredBlockStorage& blockStorage, const IBlock& block );

inline void transformGlobalToBlockLocal( CellSet& local, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellSet& global );
inline void transformGlobalToBlockLocal( CellSet& cells, const StructuredBlockStorage& blockStorage, const IBlock& block );

inline void transformBlockLocalToGlobal( CellSet& global, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellSet& local );
inline void transformBlockLocalToGlobal( CellSet& cells,  const StructuredBlockStorage& blockStorage, const IBlock& block );



// implementations

/// global cells are transformed to the block local cell space and added to vector 'local' via calling push_back !

inline void transformGlobalToBlockLocal( CellVector& local, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellVector& global ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellVector::const_iterator it = global.begin(); it != global.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformGlobalToBlockLocalCell( cell, block );
      local.push_back( cell );
   }
}



/// cells in vector 'cells' are transformed in place (from the global to the block local cell space)

inline void transformGlobalToBlockLocal( CellVector& cells, const StructuredBlockStorage& blockStorage, const IBlock& block ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellVector::iterator cell = cells.begin(); cell != cells.end(); ++cell )
      blockStorage.transformGlobalToBlockLocalCell( *cell, block );
}



/// block local cells are transformed to the global cell space and added to vector 'global' via calling push_back !

inline void transformBlockLocalToGlobal( CellVector& global, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellVector& local ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellVector::const_iterator it = local.begin(); it != local.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformBlockLocalToGlobalCell( cell, block );
      global.push_back( cell );
   }
}



/// cells in vector 'cells' are transformed in place (from the block local to the global cell space)

inline void transformBlockLocalToGlobal( CellVector& cells, const StructuredBlockStorage& blockStorage, const IBlock& block ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellVector::iterator cell = cells.begin(); cell != cells.end(); ++cell )
      blockStorage.transformBlockLocalToGlobalCell( *cell, block );
}



/// global cells are transformed to the block local cell space and added to set 'local' via calling insert ! [-> O(N*logN)]

inline void transformGlobalToBlockLocal( CellSet& local, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellSet& global ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellSet::const_iterator it = global.begin(); it != global.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformGlobalToBlockLocalCell( cell, block );
      local.insert( cell );
   }
}



/// cells in set 'cells' are transformed in place (from the global to the block local cell space) [cells are possibly reordered -> O(N*logN)]

inline void transformGlobalToBlockLocal( CellSet& cells, const StructuredBlockStorage& blockStorage, const IBlock& block ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   CellSet localCells;

   for( CellSet::const_iterator it = cells.begin(); it != cells.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformGlobalToBlockLocalCell( cell, block );
      localCells.insert( localCells.end(), cell );
   }

   cells.swap( localCells );
}



/// block local cells are transformed to the global cell space and added to set 'global' via calling insert ! [-> O(N*logN)]

inline void transformBlockLocalToGlobal( CellSet& global, const StructuredBlockStorage& blockStorage, const IBlock& block, const CellSet& local ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   for( CellSet::const_iterator it = local.begin(); it != local.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformBlockLocalToGlobalCell( cell, block );
      global.insert( cell );
   }
}



/// cells in set 'cells' are transformed in place (from the block local to the global cell space) [cells are possibly reordered -> O(N*logN)]

inline void transformBlockLocalToGlobal( CellSet& cells, const StructuredBlockStorage& blockStorage, const IBlock& block ) {

   WALBERLA_ASSERT_EQUAL( &( blockStorage.getBlockStorage() ), &( block.getBlockStorage() ) );

   CellSet globalCells;

   for( CellSet::const_iterator it = cells.begin(); it != cells.end(); ++it ) {
      Cell cell = *it;
      blockStorage.transformBlockLocalToGlobalCell( cell, block );
      globalCells.insert( globalCells.end(), cell );
   }

   cells.swap( globalCells );
}



} // namespace domain_decomposition
} // namespace walberla
