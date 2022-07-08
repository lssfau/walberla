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
//! \file GetOredNeighborhood.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Combines the flags of all neighboring cells using bitwise OR.
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/timing/TimingPool.h"

#include "field/FlagField.h"

#include <type_traits>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Combines the flags of all neighboring cells using bitwise OR.
 * Flags are read from flagField and stored in neighborhoodFlagField. Every cell contains the bitwise OR of the
 * neighboring flags in flagField.
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T >
void getOredNeighborhood(const FlagField_T* flagField, FlagField_T* neighborhoodFlagField)
{
   WALBERLA_ASSERT_GREATER_EQUAL(flagField->nrOfGhostLayers(), 2);
   WALBERLA_ASSERT_EQUAL(neighborhoodFlagField->xyzSize(), flagField->xyzSize());
   WALBERLA_ASSERT_EQUAL(neighborhoodFlagField->xyzAllocSize(), flagField->xyzAllocSize());

   // REMARK: here is the reason why the flag field MUST have two ghost layers;
   // the "OredNeighborhood" of the first ghost layer is determined, such that the first ghost layer's neighbors (i.e.,
   // the second ghost layer) must be available
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField, uint_c(1), {
      const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);
      const typename FlagField_T::Ptr neighborhoodFlagFieldPtr(*neighborhoodFlagField, x, y, z);

      *neighborhoodFlagFieldPtr = field::getOredNeighborhood< Stencil_T >(flagFieldPtr);
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOSTLAYER_XYZ
}

} // namespace free_surface
} // namespace walberla
