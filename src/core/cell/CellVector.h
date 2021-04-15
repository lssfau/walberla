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
//! \file CellVector.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Cell.h"
#include "CellInterval.h"

#include "core/DataTypes.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <vector>


namespace walberla {
namespace cell {



/// A vector of cells

class CellVector
{
public:

   /*! \name Standard container typedefs */
   //@{
   using iterator = std::vector<Cell>::iterator;
   using const_iterator = std::vector<Cell>::const_iterator;
   using size_type = std::vector<Cell>::size_type;
   using reference = std::vector<Cell>::reference;
   using const_reference = std::vector<Cell>::const_reference;
   using difference_type = std::vector<Cell>::difference_type;
   using value_type = Cell;
   //@}

   /*! \name Constructors */
   //@{
   CellVector() = default;
   CellVector(size_type n, const Cell & value = Cell()) : cells_(n, value) {}
   template <class InputIterator>
      CellVector(InputIterator first, InputIterator last) : cells_(first, last) { }
   //@}

   /*! \name Vector functions
    *  \brief see documentation of std::vector */
   //@{
   iterator       begin()        { return cells_.begin(); }
   const_iterator begin()  const { return cells_.begin(); }
   const_iterator cbegin() const { return cells_.begin(); }

   iterator       end()        { return cells_.end(); }
   const_iterator end()  const { return cells_.end(); }
   const_iterator cend() const { return cells_.end(); }

   reference       operator[](size_type n)       { return cells_[n]; }
   const_reference operator[](size_type n) const { return cells_[n]; }

   reference       at(size_type n)       { return cells_.at(n); }
   const_reference at(size_type n) const { return cells_.at(n); }

   reference       front()       { return cells_.front(); }
   const_reference front() const { return cells_.front(); }

   reference       back()       { return cells_.back(); }
   const_reference back() const { return cells_.back(); }

   void reserve(size_type n) { cells_.reserve(n); }
   void resize(size_type n, const Cell & value = Cell()) { cells_.resize(n, value); }

   size_type size() const { return cells_.size(); }
   bool empty() const { return cells_.empty(); }

   template <class InputIterator>
   void assign(InputIterator first, InputIterator last) { cells_.assign(first, last); }

   void push_back(const Cell & x) { cells_.push_back(x); }

   void pop_back() { cells_.pop_back(); }

   iterator erase(iterator position) { return cells_.erase(position);}
   iterator erase(iterator first, iterator last) { return cells_.erase(first, last);}

   iterator insert(iterator position, const Cell & cell) { return cells_.insert(position, cell);}
   void     insert(iterator position, size_type n, const Cell & cell) { cells_.insert(position, n, cell); }
   template <class InputIterator>
      void  insert(iterator position, InputIterator first, InputIterator last) { cells_.insert(position, first, last); }

   void swap(CellVector & rhs) { cells_.swap(rhs.cells_); }
   void clear() { cells_.clear(); }
   //@}

   /*! \name Utility functions */
   //@{
   inline void push_back( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ); ///< push_back( Cell(x,y,z) )
   inline void push_back( const uint_t     x, const uint_t     y, const uint_t     z ); ///< push_back( Cell(x,y,z) )

   difference_type removeDuplicates();

   CellInterval boundingBox() const;

   inline bool contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline bool contains( const uint_t     x, const uint_t     y, const uint_t     z ) const;

   inline void sort() { std::sort(cells_.begin(), cells_.end()); }
   //@}

   friend inline bool operator==( const CellVector& lhs, const CellVector& rhs ) { return lhs.cells_ == rhs.cells_; } ///< compares the content of two CellVectors
   friend inline bool operator!=( const CellVector& lhs, const CellVector& rhs ) { return lhs.cells_ != rhs.cells_; } ///< compares the content of two CellVectors

protected:

   std::vector<Cell> cells_;  ///< The vector the contained cells are stored in

}; // class CellVector


std::ostream & operator<<( std::ostream & os, const CellVector & cells );


inline void CellVector::push_back( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) {

   cells_.emplace_back(x, y, z );
}



inline void CellVector::push_back( const uint_t x, const uint_t y, const uint_t z ) {

   cells_.emplace_back(x, y, z );
}



//**********************************************************************************************************************
/*! \brief Return if cell (x,y,z) is contained in cell vector
 *
 * Complexity is O(N), where N == this->size(). If you need a data structure for storing cells that
 * provides an contains method that runs in O(logN) use a CellSet, not a CellVector.
 *
 * \return true, if cell is contained
 */
//**********************************************************************************************************************
inline bool CellVector::contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const {

   return std::find( cells_.begin(), cells_.end(), Cell(x,y,z) ) != cells_.end();
}



/// for documentation see "contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )"
inline bool CellVector::contains( const uint_t x, const uint_t y, const uint_t z ) const {

   return std::find( cells_.begin(), cells_.end(), Cell(x,y,z) ) != cells_.end();
}



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
          typename G >    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const CellVector & cellVector )
{
   mpi::sendContainer(buf, cellVector);
   return buf;
}

template< typename T >    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, CellVector & cellVector )
{
   mpi::recvContainer(buf, cellVector);
   return buf;
}



} // namespace cell

using cell::CellVector;

} // namespace walberla
