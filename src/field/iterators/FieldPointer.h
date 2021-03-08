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
//! \file FieldPointer.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/cell/Cell.h"
#include "stencil/Directions.h"


namespace walberla {
namespace field {

   //*******************************************************************************************************************
   /*!
    * \brief Class representing on Field coordinate
    *
    * \ingroup field
    *
    * What is this class good for?
    *
    * Suppose there is a PDF field, and a function that calculates the density of a cell
    * \code
    *    void calcDensity( const GhostLayerField<real_t,19>::base_iterator & it ) {
    *       // calculates density using it.neighbor() or it.getF()
    *    }
    * \endcode
    *
    * This function cannot be called when only a x,y,z coordinate is available.
    * A FieldPointer is a class that behaves like a FieldIterator but can be constructed from
    * x,y,z,f coordinates. So the code can be rewritten as:
    * \code
    *    template< typename Iter>
    *    void calcDensity( const Iter & it );
    *
    *    // and calculate the density of cell 3,2,4 using a FieldPointer
    *    calcDensity( GhostLayerField<real_t,19>::Ptr( 3,2,4) );
    * \endcode
    *
    *
    * Do not use this directly use the typedef Field::Ptr
    *
    */
   //*******************************************************************************************************************
   template<typename Field_T, typename FieldMember, typename Value_T>
   class FieldPointer
   {
   public:
      using value_type = Value_T;
      static const uint_t F_SIZE = Field_T::F_SIZE;

      FieldPointer( FieldMember & field, cell_idx_t _x, cell_idx_t _y, cell_idx_t _z, cell_idx_t _f = 0 )
         : x_(_x), y_(_y), z_(_z), f_(_f), field_( field )
      {}

      void set( cell_idx_t _x, cell_idx_t _y, cell_idx_t _z ) {
         x_ = _x;
         y_ = _y;
         z_ = _z;
      }

      inline value_type & operator*()  const { return field_(x_,y_,z_,f_); }
      inline value_type & operator->() const { return field_(x_,y_,z_,f_); }

      inline value_type & getF( cell_idx_t cf ) const { return field_(x_,y_,z_,cf); }
      inline value_type & getF( uint_t cf )     const { return field_(x_,y_,z_,cf); }

      inline value_type & operator[]( cell_idx_t cf ) const { return field_(x_,y_,z_,cf); }
      inline value_type & operator[]( uint_t cf )     const { return field_(x_,y_,z_,cf); }

      inline value_type & neighbor( stencil::Direction d, cell_idx_t cf = 0 ) const
      {
         return field_( x_ + stencil::cx[d],
                        y_ + stencil::cy[d],
                        z_ + stencil::cz[d],
                        f_ + cf );
      }

      inline value_type & neighbor( stencil::Direction d, uint_t cf ) const {
         return neighbor( d, cell_idx_c(cf) );
      }

      inline value_type & neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, cell_idx_t cf = 0 ) const {
         return field_( x_ + cx, y_ + cy,  z_ + cz, f_ + cf );
      }

      inline value_type & neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, uint_t cf ) const {
         return neighbor( cx,cy,cz, cell_idx_c(cf) );
      }

      inline cell_idx_t & x()       { return x_; }
      inline cell_idx_t   x() const { return x_; }
      inline cell_idx_t & y()       { return y_; }
      inline cell_idx_t   y() const { return y_; }
      inline cell_idx_t & z()       { return z_; }
      inline cell_idx_t   z() const { return z_; }
      inline cell_idx_t & f()       { return f_; }
      inline cell_idx_t   f() const { return f_; }

      inline Cell cell() const { return Cell( x_, y_, z_ ); }

      const FieldMember * getField() const { return &field_; }

   protected:

      cell_idx_t x_;
      cell_idx_t y_;
      cell_idx_t z_;
      cell_idx_t f_;
      FieldMember & field_;
   };



} // namespace field
} // namespace walberla
