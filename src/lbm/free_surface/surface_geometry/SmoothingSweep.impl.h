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
//! \file SmoothingSweep.impl.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Smooth fill levels (used for finite difference curvature computation).
//
//======================================================================================================================

#include "core/debug/CheckFunctions.h"
#include "core/math/Utility.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include <algorithm>
#include <cmath>

#include "ContactAngle.h"
#include "SmoothingSweep.h"
#include "Utility.h"

namespace walberla
{
namespace free_surface
{
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void SmoothingSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(IBlock* const block)
{
   // get fields
   ScalarField_T* const smoothFillField = block->getData< ScalarField_T >(smoothFillFieldID_);
   const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);
   const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);

   // get flags
   const flag_t liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);
   const flag_t obstacleFlagMask           = flagField->getMask(obstacleFlagIDSet_);

   // const KernelK8< Stencil_T > smoothingKernel(real_c(2.0));

   const uint_t kernelSize = smoothingKernel_.getStencilSize();

   WALBERLA_CHECK_GREATER_EQUAL(
      smoothFillField->nrOfGhostLayers(), kernelSize,
      "Support radius of smoothing kernel results in a smoothing stencil size that exceeds the ghost layers.");

   // including ghost layers is not necessary, even if solid cells are located in the (outermost global) ghost layer,
   // because fill level in obstacle cells is set by ObstacleFillLevelSweep
   WALBERLA_FOR_ALL_CELLS(smoothFillFieldIt, smoothFillField, fillFieldIt, fillField, flagFieldIt, flagField, {
      // IMPORTANT REMARK: do not restrict this algorithm to interface cells and their neighbors; when the normals are
      // computed in the neighborhood of interface cells, the second neighbors of interface cells are also required

      // mollify fill level in interface, liquid, and gas according to equation (9) in Williams et al.
      if (isPartOfMaskSet(flagFieldIt, liquidInterfaceGasFlagMask))
      {
         real_t normalizationConstant = real_c(0);
         real_t smoothedFillLevel     = real_c(0);
         if constexpr (Stencil_T::D == uint_t(2))
         {
            for (int j = -int_c(kernelSize); j <= int_c(kernelSize); ++j)
            {
               for (int i = -int_c(kernelSize); i <= int_c(kernelSize); ++i)
               {
                  const Vector3< real_t > dirVector(real_c(i), real_c(j), real_c(0));

                  if (isPartOfMaskSet(flagFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(0)),
                                      obstacleFlagMask))
                  {
                     if (includeObstacleNeighbors_)
                     {
                        // in solid cells, use values from smoothed fill field (instead of regular fill field) that have
                        // been set by ObstacleFillLevelSweep
                        smoothedFillLevel += smoothingKernel_.kernelFunction(dirVector) *
                                             smoothFillFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(0));

                        if (dirVector.length() < smoothingKernel_.getSupportRadius())
                        {
                           normalizationConstant += smoothingKernel_.kernelFunction(dirVector);
                        }
                     } // else: do not include this direction in smoothing
                  }
                  else
                  {
                     smoothedFillLevel += smoothingKernel_.kernelFunction(dirVector) *
                                          fillFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(0));

                     if (dirVector.length() < smoothingKernel_.getSupportRadius())
                     {
                        normalizationConstant += smoothingKernel_.kernelFunction(dirVector);
                     }
                  }
               }
            }
         }
         else
         {
            if constexpr (Stencil_T::D == uint_t(3))
            {
               for (int k = -int_c(kernelSize); k <= int_c(kernelSize); ++k)
               {
                  for (int j = -int_c(kernelSize); j <= int_c(kernelSize); ++j)
                  {
                     for (int i = -int_c(kernelSize); i <= int_c(kernelSize); ++i)
                     {
                        const Vector3< real_t > dirVector(real_c(i), real_c(j), real_c(k));

                        if (isPartOfMaskSet(flagFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(k)),
                                            obstacleFlagMask))
                        {
                           if (includeObstacleNeighbors_)
                           {
                              // in solid cells, use values from smoothed fill field (instead of regular fill field)
                              // that have been set by ObstacleFillLevelSweep
                              smoothedFillLevel +=
                                 smoothingKernel_.kernelFunction(dirVector) *
                                 smoothFillFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(k));

                              if (dirVector.length() < smoothingKernel_.getSupportRadius())
                              {
                                 normalizationConstant += smoothingKernel_.kernelFunction(dirVector);
                              }
                           } // else: do not include this direction in smoothing
                        }
                        else
                        {
                           smoothedFillLevel += smoothingKernel_.kernelFunction(dirVector) *
                                                fillFieldIt.neighbor(cell_idx_c(i), cell_idx_c(j), cell_idx_c(k));

                           if (dirVector.length() < smoothingKernel_.getSupportRadius())
                           {
                              normalizationConstant += smoothingKernel_.kernelFunction(dirVector);
                           }
                        }
                     }
                  }
               }
            }
         }

         smoothedFillLevel /= normalizationConstant;
         *smoothFillFieldIt = smoothedFillLevel;
      }
   }) // WALBERLA_FOR_ALL_CELLS
}
} // namespace free_surface
} // namespace walberla
