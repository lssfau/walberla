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
//! \file Cell.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <algorithm>
#include <iterator>
#include <iostream>


namespace walberla {
namespace cell {



/*******************************************************************************************************************//**
 * \brief   A representation of a Cell's coordinates (in 3D)
 **********************************************************************************************************************/
class Cell
{
public:

   /*! \name Constructors */
   //@{
   Cell() {}
   inline Cell( const cell_idx_t _x, const cell_idx_t _y, const cell_idx_t _z ) { cell[0] = _x; cell[1] = _y; cell[2] = _z; }
 //inline Cell( const int        _x, const int        _y, const int        _z );
   inline Cell( const uint_t     _x, const uint_t     _y, const uint_t     _z );
   //@}

   /*! \name Arithmetic operators */
   //@{
   inline Cell operator+( const Cell & rhs ) const;
   inline Cell operator-( const Cell & rhs ) const;

   inline Cell & operator+=( const Cell & rhs );
   inline Cell & operator-=( const Cell & rhs );

   inline Cell operator+() const;
   inline Cell operator-() const;
   //@}

   /*! \name Comparison operators */
   //@{
   bool operator< ( const Cell & rhs ) const;
   bool operator==( const Cell & rhs ) const;
   bool operator!=( const Cell & rhs ) const { return !operator==( rhs ); }
   //@}

   /*! \name Access operators */
   //@{
   cell_idx_t   operator[]( size_t idx ) const;
   cell_idx_t & operator[]( size_t idx );

   cell_idx_t   x() const { return cell[0]; }
   cell_idx_t & x()       { return cell[0]; }
   cell_idx_t   y() const { return cell[1]; }
   cell_idx_t & y()       { return cell[1]; }
   cell_idx_t   z() const { return cell[2]; }
   cell_idx_t & z()       { return cell[2]; }
   //@}

   bool positiveIndicesOnly() const { return x() >= cell_idx_c(0) && y() >= cell_idx_c(0) && z() >= cell_idx_c(0); }

private:

   cell_idx_t cell[3]; ///< Array of the cells coordinates. cell == {x, y, z}.
};



/*! \name Stream operators */
//@{
std::ostream & operator<<( std::ostream & os, const Cell & cell );
std::istream & operator>>( std::istream & is,       Cell & cell );
//@}



// inline Cell::Cell( const int _x, const int _y, const int _z ) {
//
//    x() = cell_idx_c( _x ); y() = cell_idx_c( _y ); z() = cell_idx_c( _z );
// }



inline Cell::Cell( const uint_t _x, const uint_t _y, const uint_t _z )
{
   cell[0] = cell_idx_c( _x );
   cell[1] = cell_idx_c( _y );
   cell[2] = cell_idx_c( _z );
}



/*******************************************************************************************************************//**
 * \brief   Less-than comparison operator for Cells.
 *
 * Compares a cell's coordinates lexicographically (first x, than eventualy y and (if necessary) finally z).
 *
 * \param [in] rhs  the cell compared to *this.
 *
 * \return  \code std::lexicographical_compare(this->cell, this->cell + 3, rhs.cell, rhs.cell + 3). \endcode
 **********************************************************************************************************************/
inline bool Cell::operator<( const Cell & rhs ) const
{
   return std::lexicographical_compare( std::reverse_iterator<const cell_idx_t*>( this->cell + 3 ), std::reverse_iterator<const cell_idx_t*>( this->cell ),
                                        std::reverse_iterator<const cell_idx_t*>( rhs.cell + 3 ),   std::reverse_iterator<const cell_idx_t*>( rhs.cell ) );
}



 /******************************************************************************************************************//**
 * \brief   Equal comparison operator for Cells.
 *
 * Compares a cell's coordinates for equality (first x, then eventualy y and (if necessary) finally z).
 *
 * \param [in] rhs  The cell compared to *this.
 *
 * \return  \code (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z) \endcode
 **********************************************************************************************************************/
inline bool Cell::operator==( const Cell & rhs ) const
{
   return ( this->x() == rhs.x() ) && ( this->y() == rhs.y() ) && ( this->z() == rhs.z() );
}



/*******************************************************************************************************************//**
 * \brief   Operator providing read-only element access of a const Cell.
 *
 * \param [in] idx  Zero-based index of the cell's coordinate component.
 *
 * \return  The idx-th coordinate component. This is equal to this->cell[i].
 **********************************************************************************************************************/
inline cell_idx_t Cell::operator[]( size_t idx ) const
{
   WALBERLA_ASSERT_LESS( idx, 3, "Index 'idx' = " << idx << " out of bounds! Cell: " << *this );
   return cell[idx];
}

/*******************************************************************************************************************//**
 * \brief   Operator for component-wise addition of two cells.
 *
 * \param [in] rhs  The cell added to this.
 *
 * \return  a Cell which components are the sum of this and rhs components
 **********************************************************************************************************************/
Cell Cell::operator+( const Cell & rhs ) const
{
   Cell result;
   result.cell[0] = cell[0] + rhs.cell[0];
   result.cell[1] = cell[1] + rhs.cell[1];
   result.cell[2] = cell[2] + rhs.cell[2];
   return result;
}

/*******************************************************************************************************************//**
 * \brief   Operator for component-wise subtraction of two cells.
 *
 * \param [in] rhs  The cell subtracted from this.
 *
 * \return  a Cell which components are the difference of this and rhs components
 **********************************************************************************************************************/
Cell Cell::operator-( const Cell & rhs ) const
{
   Cell result;
   result.cell[0] = cell[0] - rhs.cell[0];
   result.cell[1] = cell[1] - rhs.cell[1];
   result.cell[2] = cell[2] - rhs.cell[2];
   return result;
}

/*******************************************************************************************************************//**
 * \brief   Additive compound assignment operator.
 *
 * \param [in] rhs  The cell added to this.
 *
 * \return  Reference to this.
 **********************************************************************************************************************/
Cell & Cell::operator+=( const Cell & rhs )
{
   cell[0] += rhs.cell[0];
   cell[1] += rhs.cell[1];
   cell[2] += rhs.cell[2];

   return *this;
}

/*******************************************************************************************************************//**
 * \brief   Subtractive compound assignment operator.
 *
 * \param [in] rhs  The cell subtracted from this.
 *
 * \return  Reference to this.
 **********************************************************************************************************************/
Cell & Cell::operator-=( const Cell & rhs )
{
   cell[0] -= rhs.cell[0];
   cell[1] -= rhs.cell[1];
   cell[2] -= rhs.cell[2];

   return *this;
}

/*******************************************************************************************************************//**
 * \brief   Unary plus operator.
 *
 * \return  *this unmodified.
 **********************************************************************************************************************/
inline Cell Cell::operator+() const
{
   return *this;
}

/*******************************************************************************************************************//**
 * \brief   Unary negation operator.
 *
 * \return  Cell with negated components.
 **********************************************************************************************************************/
inline Cell Cell::operator-() const
{
   return Cell( -x(), -y(), -z() );
}



 /******************************************************************************************************************//**
 * \brief   Operator providing element access of a Cell.
 *
 * \param [in] idx  Zero-based index of the cell's coordinate component.
 *
 * \return  The idx-th coordinate component. This is equal to this->cell[i].
 **********************************************************************************************************************/
inline cell_idx_t & Cell::operator[]( size_t idx )
{
   WALBERLA_ASSERT_LESS( idx, 3, "Index 'idx' = " << idx << " out of bounds! Cell: " << *this );
   return cell[idx];
}



/*******************************************************************************************************************//**
 * \brief   Stream output operator for a Cell object.
 *
 * Serializes a Cell like (x,y,z). Example: (7,36,1211).
 *
 * \param [in,out] os  The output stream.
 * \param [in]   cell  The cell to be serialized.
 *
 * \return  A reference to the modified output stream.
 **********************************************************************************************************************/
inline std::ostream & operator<<( std::ostream & os, const Cell & cell )
{
   os << "(" << cell.x() << "," << cell.y() << "," << cell.z() << ")";
   return os;
}

/*******************************************************************************************************************//**
 * \brief   Stream input operator for a Cell object.
 *
 * Deserializes a Cell like (x,y,z). Example: (7,36,1211). Whitespaces between commas or brackets
 * are ignored.
 *
 * \param [in,out] is  The input stream.
 * \param [out]  cell  The result cell.
 *
 * \return  A reference to the modified input stream.
 **********************************************************************************************************************/
inline std::istream & operator>>( std::istream & is, Cell & cell )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1, comma2;
   cell_idx_t x(0), y(0), z(0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the vector
   if( !(is >> bracket1 >> x >> comma1 >> y >> comma2 >> z >> bracket2) ||
       bracket1 != '(' || comma1 != ',' || comma2 != ',' || bracket2 != ')' )
   {
         is.clear();
         is.seekg( pos );
         is.setstate( std::istream::failbit );
         is.flags( oldFlags );
         return is;
   }

   // Transferring the input to the vector values
   cell.x() = x; cell.y() = y; cell.z() = z;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}


/*******************************************************************************************************************//**
 * \brief   Provides a hash value for a Cell based on its coordinates.
 *
 * \param [in]   cell  The cell to be hashed.
 *
 * \return  a hopefully unique hash.
 **********************************************************************************************************************/
inline std::size_t hash_value( const Cell & cell )
{
  std::size_t seed = 0;
  std::hash<cell_idx_t> hasher;

  seed ^= hasher(cell.x()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  seed ^= hasher(cell.y()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  seed ^= hasher(cell.z()) + 0x9e3779b9 + (seed<<6) + (seed>>2);

  return seed;
}



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
          typename G >   // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Cell & cell )
{
   return buf << cell.x() << cell.y() << cell.z();
}

template< typename T >   // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Cell & cell )
{
   return buf >> cell.x() >> cell.y() >> cell.z();
}


} // namespace cell

using cell::Cell;

namespace mpi
{
   template<>
   struct BufferSizeTrait< Cell> {
      static const bool constantSize = true;
      static const uint_t size = 3 * sizeof(cell_idx_t) + BUFFER_DEBUG_OVERHEAD;
   };
}


} // namespace walberla

namespace std
{
   template<>
   struct hash< walberla::Cell >
   {
      std::size_t operator()( walberla::Cell const & cell ) const noexcept
      {
         return walberla::cell::hash_value( cell );
      }
   };
} // namespace std
