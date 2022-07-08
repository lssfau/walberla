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
//! \file ObstacleFillLevelSweep.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Reflect fill levels into obstacle cells (for finite difference curvature computation).
//
//======================================================================================================================

#include "core/debug/CheckFunctions.h"
#include "core/math/Utility.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include <algorithm>
#include <cmath>

namespace walberla
{
namespace free_surface
{
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ObstacleFillLevelSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(IBlock* const block)
{
   // get fields
   const ScalarField_T* const fillFieldSrc        = block->getData< const ScalarField_T >(fillFieldSrcID_);
   ScalarField_T* const fillFieldDst              = block->getData< ScalarField_T >(fillFieldDstID_);
   const FlagField_T* const flagField             = block->getData< const FlagField_T >(flagFieldID_);
   const VectorField_T* const obstacleNormalField = block->getData< const VectorField_T >(obstacleNormalFieldID_);

   // get flags
   const flag_t liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);
   const flag_t obstacleFlagMask           = flagField->getMask(obstacleFlagIDSet_);

   // equation (4.22) in dissertation of S. Bogner, 2017 (section 4.4.2.1); include ghost layer because solid cells
   // might be located in the (outermost global) ghost layer
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(fillFieldDst, uint_c(1), {
      // IMPORTANT REMARK: do not restrict this algorithm to obstacle cells that are direct neighbors of interface
      // cells; the succeeding SmoothingSweep uses the values computed here and must be executed for an at least
      // two-cell neighborhood of interface cells
      const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);
      const typename ScalarField_T::ConstPtr fillFieldSrcPtr(*fillFieldSrc, x, y, z);
      const typename ScalarField_T::Ptr fillFieldDstPtr(*fillFieldDst, x, y, z);

      if (isPartOfMaskSet(flagFieldPtr, obstacleFlagMask) &&
          isFlagInNeighborhood< Stencil_T >(flagFieldPtr, liquidInterfaceGasFlagMask))
      {
         WALBERLA_CHECK_GREATER(obstacleNormalField->get(x, y, z).length(), real_c(0),
                                "An obstacleNormal of an obstacle cell was found to be zero in obstacleNormalSweep. "
                                "This is not plausible.");

         real_t sum       = real_c(0);
         real_t weightSum = real_c(0);
         for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
         {
            if (isPartOfMaskSet(flagFieldPtr.neighbor(*dir), liquidInterfaceGasFlagMask))
            {
               const Vector3< real_t > dirVector =
                  Vector3< real_t >(real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz())).getNormalized();

               const real_t weight = std::abs(obstacleNormalField->get(x, y, z) * dirVector);

               sum += weight * fillFieldSrcPtr.neighbor(*dir);

               weightSum += weight;
            }
         }

         *fillFieldDstPtr = weightSum > real_c(0) ? sum / weightSum : real_c(0);
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
} // namespace free_surface
} // namespace walberla
