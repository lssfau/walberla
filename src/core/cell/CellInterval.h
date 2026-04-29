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
//! \file CellInterval.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include "Cell.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <algorithm>
#include <iterator>


namespace walberla {
namespace cell {



class CellSet;
class CellVector;



/// CellInterval represents a cell bounding box

class CellIntervalIterator;

class CellInterval {

public:
   using const_iterator = CellIntervalIterator;

   WALBERLA_HOST_DEVICE CellInterval() : min_( cell_idx_c(0), cell_idx_c(0), cell_idx_c(0) ), max_( cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(-1) ) {}
   WALBERLA_HOST_DEVICE CellInterval( const Cell& _min, const Cell& _max ) : min_( _min ), max_( _max ) {}
   WALBERLA_HOST_DEVICE CellInterval( const cell_idx_t _xMin, const cell_idx_t _yMin, const cell_idx_t _zMin,
                 const cell_idx_t _xMax, const cell_idx_t _yMax, const cell_idx_t _zMax ) :
                 min_( _xMin, _yMin, _zMin ), max_( _xMax, _yMax, _zMax ) {}

   WALBERLA_HOST_DEVICE bool operator==( const CellInterval& rhs ) const { return ( empty() && rhs.empty() ) || ( min_ == rhs.min_ && max_ == rhs.max_ ); }
   WALBERLA_HOST_DEVICE bool operator!=( const CellInterval& rhs ) const { return !operator==( rhs ); }

   WALBERLA_HOST_DEVICE inline Cell& min() { return min_; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& xMin() { return min_[0]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& yMin() { return min_[1]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& zMin() { return min_[2]; }

   WALBERLA_HOST_DEVICE inline Cell& max() { return max_; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& xMax() { return max_[0]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& yMax() { return max_[1]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t& zMax() { return max_[2]; }

   WALBERLA_HOST_DEVICE inline const Cell& min() const { return min_; }
   WALBERLA_HOST_DEVICE inline cell_idx_t xMin() const { return min_[0]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t yMin() const { return min_[1]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t zMin() const { return min_[2]; }

   WALBERLA_HOST_DEVICE inline const Cell& max() const { return max_; }
   WALBERLA_HOST_DEVICE inline cell_idx_t xMax() const { return max_[0]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t yMax() const { return max_[1]; }
   WALBERLA_HOST_DEVICE inline cell_idx_t zMax() const { return max_[2]; }

   WALBERLA_HOST_DEVICE inline bool empty() const { return min_.x() > max_.x() || min_.y() > max_.y() || min_.z() > max_.z(); }
   WALBERLA_HOST_DEVICE inline bool positiveIndicesOnly() const;

   WALBERLA_HOST_DEVICE inline bool contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   WALBERLA_HOST_DEVICE inline bool contains( const uint_t     x, const uint_t     y, const uint_t     z ) const { return contains( Cell(x,y,z) ); }
   WALBERLA_HOST_DEVICE inline bool contains( const Cell& cell ) const { return contains( cell.x(), cell.y(), cell.z() ); }
   WALBERLA_HOST_DEVICE inline bool contains( const CellInterval& other ) const;

   WALBERLA_HOST_DEVICE inline bool overlaps( const CellInterval& other )      const;
                               bool overlaps( const CellSet&      cellSet )    const;
                               bool overlaps( const CellVector&   cellVector ) const;

   WALBERLA_HOST_DEVICE inline CellInterval& shift( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ); ///< diagonal shift
   WALBERLA_HOST_DEVICE inline CellInterval& shift( const uint_t     x, const uint_t     y, const uint_t     z ); ///< diagonal shift
   WALBERLA_HOST_DEVICE inline CellInterval& shift( const Cell & offset );

   WALBERLA_HOST_DEVICE inline void expand( const cell_idx_t numberOfCells );
   WALBERLA_HOST_DEVICE inline void expand( const Cell & numberOfCells );

   WALBERLA_HOST_DEVICE inline void intersect( const CellInterval& other );

   WALBERLA_HOST_DEVICE inline uint_t xSize() const { return empty() ? uint_c(0) : uint_c( max_.x() - min_.x() + cell_idx_c(1) ); }
   WALBERLA_HOST_DEVICE inline uint_t ySize() const { return empty() ? uint_c(0) : uint_c( max_.y() - min_.y() + cell_idx_c(1) ); }
   WALBERLA_HOST_DEVICE inline uint_t zSize() const { return empty() ? uint_c(0) : uint_c( max_.z() - min_.z() + cell_idx_c(1) ); }

   WALBERLA_HOST_DEVICE inline uint_t size( const uint_t i ) const { return empty() ? uint_c(0) : uint_c( max_[i] - min_[i] + cell_idx_c(1) ); }

   WALBERLA_HOST_DEVICE inline uint_t numCells() const { return xSize() * ySize() * zSize(); }

   WALBERLA_HOST_DEVICE inline const_iterator begin() const;
   WALBERLA_HOST_DEVICE inline const_iterator end()   const;

private:

   Cell min_;
   Cell max_;
};



WALBERLA_HOST_DEVICE inline bool CellInterval::positiveIndicesOnly() const {

   return cell_idx_c(0) <= min_.x() && cell_idx_c(0) <= max_.x() &&
          cell_idx_c(0) <= min_.y() && cell_idx_c(0) <= max_.y() &&
          cell_idx_c(0) <= min_.z() && cell_idx_c(0) <= max_.z();
}


class CellIntervalIterator
{
public:
   using iterator_category = std::bidirectional_iterator_tag;
   using value_type = Cell;
   using difference_type = ptrdiff_t;
   using pointer = Cell *;
   using reference = Cell &;

   WALBERLA_HOST_DEVICE CellIntervalIterator( const CellInterval & ci, const Cell & cell ) : ci_(ci), cell_( cell ) { }

   WALBERLA_HOST_DEVICE const CellIntervalIterator &  operator++()    { increment(); return *this; }
   WALBERLA_HOST_DEVICE       CellIntervalIterator    operator++(int) { CellIntervalIterator tmp = *this; increment(); return tmp; }

   WALBERLA_HOST_DEVICE const CellIntervalIterator &  operator--()    { decrement(); return *this; }
   WALBERLA_HOST_DEVICE       CellIntervalIterator    operator--(int) { CellIntervalIterator tmp = *this; decrement(); return tmp; }

   WALBERLA_HOST_DEVICE bool operator==(const CellIntervalIterator & other) const { /* TODO: Device Support */ WALBERLA_ASSERT_EQUAL( &ci_, &other.ci_ ) return cell_ == other.cell_; }
   WALBERLA_HOST_DEVICE bool operator!=(const CellIntervalIterator & other) const { /* TODO: Device Support */ WALBERLA_ASSERT_EQUAL( &ci_, &other.ci_ ) return cell_ != other.cell_; }

   WALBERLA_HOST_DEVICE       Cell   operator*()  const { /* TODO: Device Support */ WALBERLA_ASSERT( ci_.contains(cell_) ) return cell_;  }
   WALBERLA_HOST_DEVICE const Cell * operator->() const { /* TODO: Device Support */ WALBERLA_ASSERT( ci_.contains(cell_) ) return &cell_; }

private:
   WALBERLA_HOST_DEVICE inline void increment();
   WALBERLA_HOST_DEVICE inline void decrement();

   const CellInterval & ci_;
   Cell cell_;
};

WALBERLA_HOST_DEVICE void CellIntervalIterator::increment()
{
   if( ++cell_.x() > ci_.xMax() )
   {
      cell_.x() = ci_.xMin();
      if( ++cell_.y() > ci_.yMax() )
      {
         cell_.y() = ci_.yMin();
         ++cell_.z();
      }
   }
}

WALBERLA_HOST_DEVICE void CellIntervalIterator::decrement()
{
   if( --cell_.x() < ci_.xMin() )
   {
      cell_.x() = ci_.xMax();
      if( --cell_.y() < ci_.yMin() )
      {
         cell_.y() = ci_.yMax();
         --cell_.z();
      }
   }
}

WALBERLA_HOST_DEVICE CellInterval::const_iterator CellInterval::begin() const
{
   if( empty() )
      return end();
   else
      return { *this, min() };
}

WALBERLA_HOST_DEVICE CellInterval::const_iterator CellInterval::end()   const { return ++CellIntervalIterator(*this, max()); }


/// Returns true only if cell (x,y,z) is contained in the cell interval

WALBERLA_HOST_DEVICE inline bool CellInterval::contains( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return ( ( min_.x() <= x ) && ( x <= max_.x() )
         && ( min_.y() <= y ) && ( y <= max_.y() )
         && ( min_.z() <= z ) && ( z <= max_.z() ) );
}



WALBERLA_HOST_DEVICE inline bool CellInterval::contains( const CellInterval& other ) const
{
   if( empty() )
      return false;
   if( other.empty() )
      return true;

   return other.min_.x() >= min_.x() && other.min_.y() >= min_.y() && other.min_.z() >= min_.z() &&
          other.max_.x() <= max_.x() && other.max_.y() <= max_.y() && other.max_.z() <= max_.z();
}



/*******************************************************************************************************************//**
 * \brief   Query if the CellIntervals overlaps with another given CellInterval.
 *
 * CellIntervals overlap if they booth include at least one common cell.
 *
 * \param   other The other CellInterval.
 *
 * \return  true if the intervals overlap, false else.
 **********************************************************************************************************************/

WALBERLA_HOST_DEVICE inline bool CellInterval::overlaps( const CellInterval& other ) const
{
   if( empty() || other.empty() )
      return false;

   // NOLINTNEXTLINE(readability-simplify-boolean-expr)
   return !(other.min_.x() > max_.x() || other.min_.y() > max_.y() || other.min_.z() > max_.z() ||
       other.max_.x() < min_.x() || other.max_.y() < min_.y() || other.max_.z() < min_.z());
}



WALBERLA_HOST_DEVICE inline CellInterval& CellInterval::shift( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) {

   min_.x() += x;
   min_.y() += y;
   min_.z() += z;

   max_.x() += x;
   max_.y() += y;
   max_.z() += z;

   return *this;
}



WALBERLA_HOST_DEVICE inline CellInterval& CellInterval::shift( const uint_t x, const uint_t y, const uint_t z ) {

   /* TODO: Device Support */ WALBERLA_ASSERT_GREATER_EQUAL( cell_idx_c(x), 0 )
   /* TODO: Device Support */ WALBERLA_ASSERT_GREATER_EQUAL( cell_idx_c(y), 0 )
   /* TODO: Device Support */ WALBERLA_ASSERT_GREATER_EQUAL( cell_idx_c(z), 0 )

   min_.x() += cell_idx_c(x);
   min_.y() += cell_idx_c(y);
   min_.z() += cell_idx_c(z);

   max_.x() += cell_idx_c(x);
   max_.y() += cell_idx_c(y);
   max_.z() += cell_idx_c(z);

   return *this;
}



WALBERLA_HOST_DEVICE inline CellInterval& CellInterval::shift( const Cell & offset )
{
   min_ += offset;
   max_ += offset;

   return *this;
}



WALBERLA_HOST_DEVICE inline void CellInterval::expand( const cell_idx_t numberOfCells )
{
   min_.x() -= numberOfCells;
   min_.y() -= numberOfCells;
   min_.z() -= numberOfCells;

   max_.x() += numberOfCells;
   max_.y() += numberOfCells;
   max_.z() += numberOfCells;
}



WALBERLA_HOST_DEVICE inline void CellInterval::expand( const Cell & numberOfCells )
{
   min_.x() -= numberOfCells.x();
   min_.y() -= numberOfCells.y();
   min_.z() -= numberOfCells.z();

   max_.x() += numberOfCells.x();
   max_.y() += numberOfCells.y();
   max_.z() += numberOfCells.z();
}



WALBERLA_HOST_DEVICE inline void CellInterval::intersect( const CellInterval& other ) 
{
   min_.x() = (xMin() > other.xMin()) ? xMin() : other.xMin();
   min_.y() = (yMin() > other.yMin()) ? yMin() : other.yMin();
   min_.z() = (zMin() > other.zMin()) ? zMin() : other.zMin();

   max_.x() = (xMax() < other.xMax()) ? xMax() : other.xMax();
   max_.y() = (yMax() < other.yMax()) ? yMax() : other.yMax();
   max_.z() = (zMax() < other.zMax()) ? zMax() : other.zMax();
}



/// output in the form of "[(3,1,5) ... (4,5,8)]"

inline std::ostream& operator<<( std::ostream& os, const CellInterval& interval )
{
   os << "[" << interval.min() << " ... " << interval.max() << "]";

   return os;
}

inline std::istream& operator>>( std::istream& is, CellInterval& interval )
{
   if( !is ) return is;

   char bracket1, bracket2, dot1, dot2, dot3;
   Cell c1(0,0,0), c2(0,0,0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the vector
   if( !(is >> bracket1 >> c1 >> dot1 >> dot2 >> dot3 >> c2 >> bracket2) ||
      bracket1 != '[' || dot1 != '.' || dot2 != '.' || dot3 != '.' || bracket2 != ']' )
   {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transferring the input to the vector values
   interval.min() = c1;
   interval.max() = c2;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
          typename G >   // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const CellInterval & cellInterval )
{
   buf.addDebugMarker( "ci" );
   return buf << cellInterval.min() << cellInterval.max();
}

template< typename T >    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, CellInterval & cellInterval )
{
   buf.readDebugMarker( "ci" );
   return buf >> cellInterval.min() >> cellInterval.max();
}


} // namespace cell

using cell::CellInterval;

namespace mpi
{
   template<>
   struct BufferSizeTrait< CellInterval> {
      static const bool constantSize = BufferSizeTrait<Cell>::constantSize;
      static const uint_t size =  2 * BufferSizeTrait<Cell>::size + BUFFER_DEBUG_OVERHEAD;
   };
}


} // namespace walberla
