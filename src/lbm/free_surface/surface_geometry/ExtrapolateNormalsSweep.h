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
//! \file ExtrapolateNormalsSweep.h
//! \ingroup surface_geometry
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Extrapolate interface normals to neighboring cells in D3Q27 direction.
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

#include <type_traits>
#include <vector>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Approximates the normals of non-interface cells using the normals of all neighboring interface cells in D3Q27
 * direction.
 * The approximation is computed by summing the weighted normals of all neighboring interface cells. The weights are
 * chosen as in the Parker-Youngs approximation.
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
class ExtrapolateNormalsSweep
{
 protected:
   using FlagUIDSet = Set< FlagUID >;

   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   ExtrapolateNormalsSweep(const BlockDataID& normalFieldID, const ConstBlockDataID& flagFieldID,
                           const FlagUID& interfaceFlagID)
      : normalFieldID_(normalFieldID), flagFieldID_(flagFieldID), interfaceFlagID_(interfaceFlagID)
   {}

   void operator()(IBlock* const block);

 private:
   BlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;

   FlagUID interfaceFlagID_;
}; // class ExtrapolateNormalsSweep

} // namespace free_surface
} // namespace walberla

#include "ExtrapolateNormalsSweep.impl.h"