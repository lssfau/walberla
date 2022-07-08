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
//! \file GetLaplacePressure.h
//! \ingroup dynamics
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute difference in density due to Laplace pressure.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/logging/Logging.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Compute difference in density due to Laplace pressure.
 **********************************************************************************************************************/
inline real_t computeDeltaRhoLaplacePressure(real_t sigma, real_t curvature, real_t maxDeltaRho = real_c(0.1))
{
   static const real_t inv_cs2  = real_c(3); // 1.0 / (cs * cs)
   const real_t laplacePressure = real_c(2) * sigma * curvature;
   real_t deltaRho              = inv_cs2 * laplacePressure;

   if (deltaRho > maxDeltaRho)
   {
      WALBERLA_LOG_WARNING("Too large density variation of " << deltaRho << " due to laplacePressure "
                                                             << laplacePressure << " with curvature " << curvature
                                                             << ". Will be limited to " << maxDeltaRho << ".\n");
      deltaRho = maxDeltaRho;
   }
   if (deltaRho < -maxDeltaRho)
   {
      WALBERLA_LOG_WARNING("Too large density variation of " << deltaRho << " due to laplacePressure "
                                                             << laplacePressure << " with curvature " << curvature
                                                             << ". Will be limited to " << -maxDeltaRho << ".\n");
      deltaRho = -maxDeltaRho;
   }

   return deltaRho;
}

} // namespace free_surface
} // namespace walberla
