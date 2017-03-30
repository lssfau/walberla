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
//! \file FastInvSqrt.h
//! \ingroup math
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Optimized implementation of 1/sqrt( double )
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace math {

   union DoubleAndLongLong { double asReal; uint64_t asInt; };
   union FloatAndLong      { float  asReal; uint32_t asInt; };

   //*******************************************************************************************************************
   /*! Computes 1.0 / sqrt(x) with limited accuracy
   *
   *  The bigger y is ( relative to 1 ) the better the approximation
   *
   *  \tparam numIter number of iterations
   *                  higher numbers make the function slower but result more accurate
   *
   */
   //*******************************************************************************************************************
   template<unsigned int numIter>
   inline double fastInvSqrt( double y )
   {
      double yhalf = ( double )0.5 * y;
      DoubleAndLongLong u;
      u.asReal = y;
      u.asInt = 0x5fe6ec85e7de30daLL - ( u.asInt >> 1 );
      y = u.asReal;
      for ( unsigned int k=0; k < numIter; ++k )
         y = y * ( 1.5 - yhalf * y * y );
      return y;
   }


   //*******************************************************************************************************************
   /*! Float version, for documentation see above
   */
   //*******************************************************************************************************************
   template<unsigned int numIter>
   inline float fastInvSqrt( float y )
   {
      float yhalf = 0.5f * y;
      FloatAndLong u;
      u.asReal = y;
      u.asInt = 0x5f3759df - ( u.asInt >> 1 );
      y = u.asReal;
      for ( unsigned int k=0; k < numIter; ++k )
         y = y * ( 1.5f - yhalf * y * y );

      return y;
   }


} // namespace math
} // namespace walberla


