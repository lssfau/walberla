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
//! \file VectorFieldAccessor.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include <type_traits>

namespace walberla {
namespace field {


   /// Provides an abstraction between Field<real_t,3> and Field<Vector3<real_t>, 1 >
   template<typename VectorField_T, typename Enable=void >
   struct VectorFieldAccessor
   {
      static_assert( VectorField_T::F_SIZE == 3, "Only valid for Fields with 3 components (F_SIZE==3)" );
      static_assert( std::is_same< typename VectorField_T::value_type, real_t >::value, "Only works for real valued fields" );

      typedef Vector3<real_t> vector_or_constRefVector;

      static vector_or_constRefVector get( const VectorField_T * f, cell_idx_t x, cell_idx_t y, cell_idx_t z )
      {
          return Vector3<real_t>( f->get(x,y,z,0), f->get(x,y,z,1), f->get(x,y,z,2) );
      }

      static void set( const VectorField_T * f, const Vector3<real_t> & v, cell_idx_t x, cell_idx_t y, cell_idx_t z )
      {
          f->get(x,y,z,0) = v[0];
          f->get(x,y,z,1) = v[1];
          f->get(x,y,z,2) = v[2];
      }
   };

   template<typename VectorField_T>
   struct VectorFieldAccessor<VectorField_T,
                              typename std::enable_if< std::is_same< typename VectorField_T::value_type,
                                                                           Vector3<real_t> >::value >::type >
   {
       typedef const Vector3<real_t> & vector_or_constRefVector;

       static vector_or_constRefVector get( const VectorField_T * f, cell_idx_t x, cell_idx_t y, cell_idx_t z )
       {
           return f->get(x,y,z);
       }

       static void set( const VectorField_T * f, const Vector3<real_t> & v, cell_idx_t x, cell_idx_t y, cell_idx_t z )
       {
           f->get(x,y,z) = v;
       }
   };





} // namespace field
} // namespace walberla


