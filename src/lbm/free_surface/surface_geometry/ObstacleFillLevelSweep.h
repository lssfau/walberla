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

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/BlockDataID.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Reflect fill levels into obstacle cells by averaging the fill levels from fluid cells with weights according to the
 * surface normal.
 *
 * See dissertation of S. Bogner, 2017 (section 4.4.2.1).
 *
 * IMPORTANT REMARK: If an obstacle is located in a non-periodic outermost ghost layer, the fill level field must have
 * two ghost layers. That is, the ObstacleFillLevelSweep computes the fill level obstacle of cells located in the
 * outermost global ghost layer. For this, the all neighboring cells' fill levels are required.
 * A single ghost layer is not sufficient, because the computed values by ObstacleFillLevelSweep (located in an
 * outermost global ghost layer) are not communicated themselves. In the example below, the values A, D, E, and H are
 * located in a global outermost ghost layer. Only directions without # shall be communicated and * marks ghost layers
 * in the directions to be communicated. In this example, only B, C, F, and G will be communicated as expected. In
 * contrast, A, D, E, and H will not be communicated.
 *
 *  Block 1      Block 2                            Block 1      Block 2
 *  ######       ######                             ######       ######
 *  # A |*       *| E #                             # A |*       *| E #
 *  # ----       -----#                             # ----       -----#
 *  # B |*       *| F #     ===> communication      # B |F       B| F #
 *  # C |*       *| G #                             # C |G       C| G #
 *  # ----       -----#                             # ----       -----#
 *  # D |*       *| H #                             # D |*       *| H #
 *  ######       ######                             ######       ######
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ObstacleFillLevelSweep
{
 protected:
   using FlagUIDSet = Set< FlagUID >;

   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   ObstacleFillLevelSweep(const BlockDataID& fillFieldDstID, const ConstBlockDataID& fillFieldSrcID,
                          const ConstBlockDataID& flagFieldID, const ConstBlockDataID& obstacleNormalFieldID,
                          const FlagUIDSet& liquidInterfaceGasFlagIDSet, const FlagUIDSet& obstacleFlagIDSet)
      : fillFieldDstID_(fillFieldDstID), fillFieldSrcID_(fillFieldSrcID), flagFieldID_(flagFieldID),
        obstacleNormalFieldID_(obstacleNormalFieldID), liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet),
        obstacleFlagIDSet_(obstacleFlagIDSet)
   {}

   void operator()(IBlock* const block);

 private:
   BlockDataID fillFieldDstID_;
   ConstBlockDataID fillFieldSrcID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID obstacleNormalFieldID_;

   FlagUIDSet liquidInterfaceGasFlagIDSet_;
   FlagUIDSet obstacleFlagIDSet_;
}; // class ObstacleFillLevelSweep

} // namespace free_surface
} // namespace walberla

#include "ObstacleFillLevelSweep.impl.h"