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
//! \file SymmetryCheck.h
//! \ingroup field
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Field.h"
#include "core/math/Utility.h"


namespace walberla {
namespace field {


   //*******************************************************************************************************************
   /*! Test for axis symmetry in a single field
   *
   * \param field         the field to check
   * \param dimension     0 for symmetry check on x axis, 1 on y axis, 2 on z axis
   * \param fCoord        which component to check
   * \param differingCell if pointer not null and field is not symmetric, this acts as a output parameter
   *                      for the first asymmetric cell found
   */
   //*******************************************************************************************************************
   template< typename Field_T >
   bool isSymmetric( const Field_T * field, uint_t dimension, uint_t fCoord = 0, Cell * differingCell = nullptr )
   {
      WALBERLA_ASSERT_LESS( dimension, 3 );

      const Cell maxCoordValues = field->xyzSize().max();
      const uint_t coordHalf    = field->size( dimension ) / 2;

      Cell iterationLimits = maxCoordValues;
      iterationLimits[dimension] = coordHalf;

      Cell currentCell;
      for( currentCell[2] = 0; currentCell[2] <= iterationLimits[2]; ++currentCell[2] )
         for( currentCell[1] = 0; currentCell[1] <= iterationLimits[1]; ++currentCell[1] )
            for( currentCell[0] = 0; currentCell[0] <= iterationLimits[0]; ++currentCell[0] )
            {
               Cell compareCell = currentCell;
               compareCell[ dimension ] = maxCoordValues[dimension] - currentCell[dimension];

               auto value1 = field->get( currentCell, fCoord );
               auto value2 = field->get( compareCell, fCoord );
               if ( ! math::equal( value1, value2 ) )
               {
                  if ( differingCell != nullptr )
                     *differingCell = currentCell;

                  return false;
               }
            }

      return true;
   }






} // namespace field
} // namespace walberla


