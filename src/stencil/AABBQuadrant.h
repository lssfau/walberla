//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURpE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file AABBQuadrant.h
//! \ingroup stencil
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Directions.h"

#include "core/math/GenericAABB.h"

namespace walberla {
namespace stencil  {


/**
* \brief Determine sin which neighboring section of an AABB a point lies
*
* If p lies within aabb returns C. Other it is determined in which of the 26 imaginary neighboring AABBs p is located. The direction of the
* from aabb to the determined neighboring aabb is returned
*
* 2D Example:
*
* -------------------------    Let aabb be the AABB marked with asterisks. The boxes surrounding aabb are the eight (in 2D) imaginary
* |       |       |       |    neighboring boxes. Let p be marked with x. The result of the function would be NE.
* |       |       |   x   |              
* |       |       |       |
* --------*********--------
* |       *       *       |
* |       *       *       |
* |       *       *       |
* --------*********--------
* |       |       |       |
* |       |       |       |
* |       |       |       |
* -------------------------
*
*
* \tparam T The scalar data type used for the AABB and the tested point
*
* \param aabb The AABB
* \param p    The tested point
*/
template< typename U >
Direction getQuadrant( const math::GenericAABB< U > & aabb, const typename math::GenericAABB< U >::vector_type & p )
{
   if( p[0] < aabb.xMin() )
   {
      if( p[1] < aabb.yMin() )
      {
         if( p[2] < aabb.zMin() )
         {
            return BSW;
         }
         else if( p[2] < aabb.zMax() )
         {
            return SW;
         }
         else
         {
            return TSW;
         }
      }
      else if( p[1] < aabb.yMax() )
      {
         if( p[2] < aabb.zMin() )
         {
            return BW;
         }
         else if( p[2] < aabb.zMax() )
         {
            return W;
         }
         else
         {
            return TW;
         }
      }
      else
      {
         if( p[2] < aabb.zMin() )
         {
            return BNW;
         }
         else if( p[2] < aabb.zMax() )
         {
            return NW;
         }
         else
         {
            return TNW;
         }
      }
   }
   else if( p[0] < aabb.xMax() )
   {
      if( p[1] < aabb.yMin() )
      {
         if( p[2] < aabb.zMin() )
         {
            return BS;
         }
         else if( p[2] < aabb.zMax() )
         {
            return S;
         }
         else
         {
            return TS;
         }
      }
      else if( p[1] < aabb.yMax() )
      {
         if( p[2] < aabb.zMin() )
         {
            return B;
         }
         else if( p[2] < aabb.zMax() )
         {
            return C;
         }
         else
         {
            return T;
         }
      }
      else
      {
         if( p[2] < aabb.zMin() )
         {
            return BN;
         }
         else if( p[2] < aabb.zMax() )
         {
            return N;
         }
         else
         {
            return TN;
         }
      }
   }
   else
   {
      if( p[1] < aabb.yMin() )
      {
         if( p[2] < aabb.zMin() )
         {
            return BSE;
         }
         else if( p[2] < aabb.zMax() )
         {
            return SE;
         }
         else
         {
            return TSE;
         }
      }
      else if( p[1] < aabb.yMax() )
      {
         if( p[2] < aabb.zMin() )
         {
            return BE;
         }
         else if( p[2] < aabb.zMax() )
         {
            return E;
         }
         else
         {
            return TE;
         }
      }
      else
      {
         if( p[2] < aabb.zMin() )
         {
            return BNE;
         }
         else if( p[2] < aabb.zMax() )
         {
            return NE;
         }
         else
         {
            return TNE;
         }
      }
   }
}

} // namespace walberla
} // namespace stencil 
