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
//! \file CellSet.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Cell.h"
#include "CellInterval.h"
#include "CellVector.h"

#include "core/Set.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"


namespace walberla {
namespace cell {



/// A set of cells

class CellSet : public Set<Cell> {

public:

   CellSet() : Set<Cell>() {}
   CellSet( const Cell& element ) : Set<Cell>( element ) {}

   CellSet( const CellVector& cells ) {
      for( CellVector::const_iterator cell = cells.begin(); cell != cells.end(); ++cell )
         Set<Cell>::insert( *cell );
   }

   using Set<Cell>::insert;

   void insert( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) { Set<Cell>::insert( Cell(x,y,z) ); }
   void insert( const uint_t     x, const uint_t     y, const uint_t     z ) { Set<Cell>::insert( Cell(x,y,z) ); }

   using Set<Cell>::contains;

   bool contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const { return Set<Cell>::contains( Cell(x,y,z) ); }
   bool contains( const uint_t     x, const uint_t     y, const uint_t     z ) const { return Set<Cell>::contains( Cell(x,y,z) ); }

   CellInterval boundingBox() const;

   void pushToCellVector( CellVector& cellVector ) const { for( auto cell = begin(); cell != end(); ++cell ) cellVector.push_back( *cell ); }

}; // class CellSet



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
   typename G >    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const CellSet & cellSet )
{
   mpi::sendAssocContainer(buf, cellSet);
   return buf;
}

template< typename T >    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, CellSet & cellSet )
{
   mpi::recvAssocContainer(buf, cellSet);
   return buf;
}



} // namespace cell

using cell::CellSet;

} // namespace walberla
