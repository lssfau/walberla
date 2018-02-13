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
//! \file IntersectionRatio.impl.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================


namespace walberla {
namespace lbm {


template <typename Body>
real_t intersectionRatioBisection( const Body & body,
                                   const Vector3<real_t> & fluidPoint,
                                   const Vector3<real_t> & direction,
                                   const real_t epsilon )
{
   WALBERLA_ASSERT( !geometry::contains( body, fluidPoint ), "fluid point: " << fluidPoint );
   WALBERLA_ASSERT(  geometry::contains( body, fluidPoint + direction ), "fluidPoint + direction: " << fluidPoint + direction );
   
   const real_t sqEpsilon         = epsilon * epsilon;
   const real_t sqDirectionLength = direction.sqrLength();
   
   real_t q( 0.5 );
   real_t qDelta( 0.5 );

   while( qDelta * qDelta * sqDirectionLength >= sqEpsilon )
   {
      qDelta *= real_t( 0.5 );
      Vector3<real_t> p = fluidPoint + q * direction;
      if( geometry::contains( body, p ) )
      {
         q -= qDelta;
      }
      else
      {
         q += qDelta;
      }
   }

   WALBERLA_ASSERT_GREATER_EQUAL( q, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( q, real_t( 1 ) );

   return q;
}


} // namespace lbm
} // namespace walberla


