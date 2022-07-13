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
//! \file ExtrapolateNormalsSweep.impl.h
//! \ingroup surface_geometry
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Extrapolate interface normals to neighboring cells in D3Q27 direction.
//
//======================================================================================================================

#include "core/math/Vector3.h"

#include "field/GhostLayerField.h"

#include "stencil/D3Q27.h"

#include "ExtrapolateNormalsSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
void ExtrapolateNormalsSweep< Stencil_T, FlagField_T, VectorField_T >::operator()(IBlock* const block)
{
   VectorField_T* const normalField   = block->getData< VectorField_T >(normalFieldID_);
   const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);

   const auto interfaceFlag = flagField->getFlag(interfaceFlagID_);

   // compute normals in interface neighboring cells, i.e., in D3Q27 direction of interface cells
   WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, {
      if (!isFlagSet(flagFieldIt, interfaceFlag) && isFlagInNeighborhood< stencil::D3Q27 >(flagFieldIt, interfaceFlag))
      {
         uint_t count     = uint_c(0);
         vector_t& normal = *normalFieldIt;
         normal.set(real_c(0), real_c(0), real_c(0));

         // approximate the normal of non-interface cells with the normal of neighboring interface cells (weights as
         // in Parker-Youngs approximation)
         for (auto i = Stencil_T::beginNoCenter(); i != Stencil_T::end(); ++i)
         {
            if (isFlagSet(flagFieldIt.neighbor(*i), interfaceFlag))
            {
               normal += real_c(stencil::gaussianMultipliers[i.toIdx()]) * normalFieldIt.neighbor(*i);
               ++count;
            }
         }

         // normalize the normal
         normal = normal.getNormalizedOrZero();
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

} // namespace free_surface
} // namespace walberla
