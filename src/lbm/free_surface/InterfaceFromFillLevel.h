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
//! \file InterfaceFromFillLevel.h
//! \ingroup free_surface
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Check whether a cell should be an interface cell according to its properties.
//
//======================================================================================================================

#pragma once

#include "core/cell/Cell.h"

#include <type_traits>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Check whether a cell is an interface cell with respect to its fill level and its direct neighborhood (liquid cells
 * can not have gas cell in their direct neighborhood; therefore, the liquid cell is forced to be an interface cell).
 **********************************************************************************************************************/
template< typename Stencil_T, typename ScalarField_T >
inline bool isInterfaceFromFillLevel(const ScalarField_T& fillField, cell_idx_t x, cell_idx_t y, cell_idx_t z)
{
   WALBERLA_ASSERT(std::is_floating_point< typename ScalarField_T::value_type >::value,
                   "Fill level field must contain floating point values.")

   real_t fillLevel = fillField.get(x, y, z);

   // this cell is regular gas cell
   if (floatIsEqual(fillLevel, real_c(0.0), real_c(1e-14))) { return false; }

   // this cell is regular interface cell
   if (fillLevel < real_c(1.0)) { return true; }

   // check this cell's direct neighborhood for gas cells (a liquid cell can not be direct neighbor to a gas cell)
   for (auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d)
   {
      // this cell has a gas cell in its direct neighborhood; it therefore must not be a liquid cell but an interface
      // cell although its fill level is 1.0
      if (fillField.get(x + d.cx(), y + d.cy(), z + d.cz()) <= real_c(0.0)) { return true; }
   }

   // this cell is a regular fluid cell
   return false;
}

template< typename Stencil_T, typename ScalarFieldIt_T >
inline bool isInterfaceFromFillLevel(const ScalarFieldIt_T& fillFieldIt)
{
   return isInterfaceFromFillLevel< Stencil_T, typename ScalarFieldIt_T::FieldType >(
      *(fillFieldIt.getField()), fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z());
}

template< typename Stencil_T, typename ScalarField >
inline bool isInterfaceFromFillLevel(const ScalarField& fillField, const Cell& cell)
{
   return isInterfaceFromFillLevel< Stencil_T >(fillField, cell.x(), cell.y(), cell.z());
}

} // namespace free_surface
} // namespace walberla
