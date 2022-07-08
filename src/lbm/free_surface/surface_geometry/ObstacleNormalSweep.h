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
//! \file ObstacleNormalSweep.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute a mean obstacle normal in interface cells near solid boundary cells.
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Compute a mean obstacle normal in interface cells near obstacle cells, and/or in obstacle cells. This reduces the
 * influence of a stair-case approximated wall in the wetting model.
 *
 * - computeInInterfaceCells: Compute the obstacle normal in interface cells. Required when using curvature computation
 *                            based on local triangulation (not with finite difference method).
 * - computeInObstacleCells: Compute the obstacle normal in obstacle cells. Required when using curvature computation
 *                           based on the finite difference method (not with local triangulation).
 *
 * Details can be found in the dissertation of S. Donath, 2011 section 6.3.5.2.
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
class ObstacleNormalSweep
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   ObstacleNormalSweep(const BlockDataID& obstacleNormalFieldID, const ConstBlockDataID& flagFieldID,
                       const FlagUID& interfaceFlagID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                       const Set< FlagUID >& obstacleFlagIDSet, bool computeInInterfaceCells,
                       bool computeInObstacleCells, bool computeInGhostLayer)
      : obstacleNormalFieldID_(obstacleNormalFieldID), flagFieldID_(flagFieldID), interfaceFlagID_(interfaceFlagID),
        liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFlagIDSet_(obstacleFlagIDSet),
        computeInInterfaceCells_(computeInInterfaceCells), computeInObstacleCells_(computeInObstacleCells),
        computeInGhostLayer_(computeInGhostLayer)
   {
      if (!computeInInterfaceCells_ && !computeInObstacleCells_)
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
            "In ObstacleNormalSweep, you specified to neither compute the obstacle normal in interface cells, nor in "
            "obstacle cells. That is, ObstacleNormalSweep will do nothing. Please check if this is what you really "
            "want.");
      }
   }

   void operator()(IBlock* const block);

 private:
   template< typename FlagFieldIt_T >
   void computeObstacleNormalInInterfaceCell(vector_t& obstacleNormal, const FlagFieldIt_T& flagFieldIt,
                                             const flag_t& validNeighborFlagMask);

   template< typename FlagFieldIt_T >
   void computeObstacleNormalInObstacleCell(vector_t& obstacleNormal, const FlagFieldIt_T& flagFieldIt,
                                            const flag_t& liquidInterfaceGasFlagMask);

   BlockDataID obstacleNormalFieldID_;
   ConstBlockDataID flagFieldID_;

   FlagUID interfaceFlagID_;
   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   Set< FlagUID > obstacleFlagIDSet_;

   bool computeInInterfaceCells_;
   bool computeInObstacleCells_;
   bool computeInGhostLayer_;
}; // class ObstacleNormalSweep

} // namespace free_surface
} // namespace walberla

#include "ObstacleNormalSweep.impl.h"