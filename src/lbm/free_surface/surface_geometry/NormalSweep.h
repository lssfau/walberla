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
//! \file NormalSweep.h
//! \ingroup surface_geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Daniela Anderl
//! \author Stefan Donath
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute interface normal.
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
 * Compute normals in interface cells by taking the derivative of the fill level field using the Parker-Youngs
 * approximation. Near boundary cells, a modified Parker-Youngs approximation is applied with the cell being shifted by
 * 0.5 away from the boundary.
 *
 * Details can be found in the Dissertation of S. Donath, page 21f.
 *
 * IMPORTANT REMARK: In this FSLBM implementation, the normal is defined to point from liquid to gas.
 *
 * More general: compute the gradient of a given scalar field on cells that are marked with a specific flag.
 **********************************************************************************************************************/

template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class NormalSweep
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using scalar_t = typename std::remove_const< typename ScalarField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   NormalSweep(const BlockDataID& normalFieldID, const ConstBlockDataID& fillFieldID,
               const ConstBlockDataID& flagFieldID, const FlagUID& interfaceFlagID,
               const Set< FlagUID >& liquidInterfaceGasFlagIDSet, const Set< FlagUID >& obstacleFlagIDSet,
               bool computeInInterfaceNeighbors, bool includeObstacleNeighbors, bool modifyNearObstacles,
               bool computeInGhostLayer)
      : normalFieldID_(normalFieldID), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet),
        obstacleFlagIDSet_(obstacleFlagIDSet), computeInInterfaceNeighbors_(computeInInterfaceNeighbors),
        includeObstacleNeighbors_(includeObstacleNeighbors), modifyNearObstacles_(modifyNearObstacles),
        computeInGhostLayer_(computeInGhostLayer)
   {}

   void operator()(IBlock* const block);

 private:
   BlockDataID normalFieldID_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;

   FlagUID interfaceFlagID_;
   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   Set< FlagUID > obstacleFlagIDSet_;

   bool computeInInterfaceNeighbors_;
   bool includeObstacleNeighbors_;
   bool modifyNearObstacles_;
   bool computeInGhostLayer_;
}; // class NormalSweep

// namespace to use these functions outside NormalSweep, e.g., in ReinitializationSweep
namespace normal_computation
{
// compute the normal using Parker-Youngs approximation (see dissertation of S. Donath, 2011, section 2.3.3.1.1)
template< typename Stencil_T, typename vector_t, typename ScalarFieldIt_T, typename FlagFieldIt_T, typename flag_t >
void computeNormal(vector_t& normal, const ScalarFieldIt_T& fillFieldIt, const FlagFieldIt_T& flagFieldIt,
                   const flag_t& validNeighborFlagMask);

// near solid boundary cells, compute a Parker-Youngs approximation around a virtual (constructed) midpoint that is
// displaced by a distance of 0.5 away from the boundary (see dissertation of S. Donath, 2011, section 6.3.5.1)
template< typename Stencil_T, typename vector_t, typename ScalarFieldIt_T, typename FlagFieldIt_T, typename flag_t >
void computeNormalNearSolidBoundary(vector_t& normal, const ScalarFieldIt_T& fillFieldIt,
                                    const FlagFieldIt_T& flagFieldIt, const flag_t& validNeighborFlagMask,
                                    const flag_t& obstacleFlagMask);
} // namespace normal_computation

} // namespace free_surface
} // namespace walberla

#include "NormalSweep.impl.h"