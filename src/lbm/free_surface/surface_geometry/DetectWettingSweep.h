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
//! \file DetectWettingSweep.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweep for detecting cells that need to be converted to interface to obtain a smooth wetting interface.
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/math/Constants.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

#include "stencil/D2Q4.h"
#include "stencil/D3Q19.h"

#include <type_traits>
#include <vector>

#include "ContactAngle.h"
#include "Utility.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Sweep for detecting interface cells that need to be created in order to obtain a smooth interface continuation in
 * case of wetting.
 *
 * See dissertation of S. Donath, 2011 section 6.3.5.3.
 **********************************************************************************************************************/
template< typename Stencil_T, typename BoundaryHandling_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
class DetectWettingSweep
{
 protected:
   using FlagUIDSet = Set< FlagUID >;

   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;

   // restrict stencil because surface continuation in corner directions is not meaningful
   using WettingStencil_T = typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q19 >::type;

 public:
   DetectWettingSweep(BlockDataID boundaryHandling, const FlagInfo< FlagField_T >& flagInfo,
                      const ConstBlockDataID& normalFieldID, const ConstBlockDataID& fillFieldID)
      : boundaryHandlingID_(boundaryHandling), flagInfo_(flagInfo), normalFieldID_(normalFieldID),
        fillFieldID_(fillFieldID)
   {}

   void operator()(IBlock* const block);

 private:
   BlockDataID boundaryHandlingID_;
   FlagInfo< FlagField_T > flagInfo_;
   ConstBlockDataID normalFieldID_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;

}; // class DetectWettingSweep

template< typename Stencil_T, typename BoundaryHandling_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
void DetectWettingSweep< Stencil_T, BoundaryHandling_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(
   IBlock* const block)
{
   // get free surface boundary handling
   BoundaryHandling_T* const boundaryHandling = block->getData< BoundaryHandling_T >(boundaryHandlingID_);

   // get fields
   const VectorField_T* const normalField = block->getData< const VectorField_T >(normalFieldID_);
   const ScalarField_T* const fillField   = block->getData< const ScalarField_T >(fillFieldID_);
   const FlagField_T* const flagField     = boundaryHandling->getFlagField();

   // get flags
   const FlagInfo< FlagField_T >& flagInfo = flagInfo_;
   using flag_t                            = typename FlagField_T::flag_t;

   const flag_t liquidInterfaceGasFlagMask = flagInfo_.liquidFlag | flagInfo_.interfaceFlag | flagInfo_.gasFlag;

   WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, fillFieldIt, fillField, {
      // skip non-interface cells
      if (!isFlagSet(flagFieldIt, flagInfo.interfaceFlag)) { continue; }

      // skip cells that have no solid cell in their neighborhood
      if (!isFlagInNeighborhood< WettingStencil_T >(flagFieldIt, flagInfo.obstacleFlagMask)) { continue; }

      // restrict maximal and minimal angle such that the surface continuation does not become too flat
      if (*fillFieldIt < real_c(0.005) || *fillFieldIt > real_c(0.995)) { continue; }

      const Vector3< real_t > interfacePointLocation = getInterfacePoint(*normalFieldIt, *fillFieldIt);

      for (auto dir = WettingStencil_T::beginNoCenter(); dir != WettingStencil_T::end(); ++dir)
      {
         const Cell neighborCell =
            Cell(flagFieldIt.x() + dir.cx(), flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
         const flag_t neighborFlag = flagField->get(neighborCell);

         // skip neighboring cells that
         // - are not liquid, gas or interface
         // - are already marked for conversion to interface due to wetting
         // IMPORTANT REMARK: It is crucial that interface cells are NOT skipped here. Since the
         // "keepInterfaceForWettingFlag" flag is cleared in all cells after performing the conversion, interface cells
         // that were converted due to wetting must still get this flag to avoid being prematurely converted back.
         if (!isPartOfMaskSet(neighborFlag, liquidInterfaceGasFlagMask) ||
             isFlagSet(neighborFlag, flagInfo.keepInterfaceForWettingFlag))
         {
            continue;
         }

         // skip neighboring cells that do not have solid cells in their neighborhood
         bool hasObstacle = false;
         for (auto dir2 = WettingStencil_T::beginNoCenter(); dir2 != WettingStencil_T::end(); ++dir2)
         {
            const Cell neighborNeighborCell =
               Cell(flagFieldIt.x() + dir.cx() + dir2.cx(), flagFieldIt.y() + dir.cy() + dir2.cy(),
                    flagFieldIt.z() + dir.cz() + dir2.cz());
            const flag_t neighborNeighborFlag = flagField->get(neighborNeighborCell);

            if (isPartOfMaskSet(neighborNeighborFlag, flagInfo.obstacleFlagMask))
            {
               hasObstacle = true;
               break; // exit dir2 loop
            }
         }
         if (!hasObstacle) { continue; }

         // check cell edges for intersection with the interface surface plane and mark the respective neighboring for
         // conversion
         // bottom south
         if ((dir.cx() == 0 && dir.cy() == -1 && dir.cz() == -1) ||
             (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == -1) || (dir.cx() == 0 && dir.cy() == -1 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // bottom north
         if ((dir.cx() == 0 && dir.cy() == 1 && dir.cz() == -1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == -1) ||
             (dir.cx() == 0 && dir.cy() == 1 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // bottom west
         if ((dir.cx() == -1 && dir.cy() == 0 && dir.cz() == -1) ||
             (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == -1) || (dir.cx() == -1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // bottom east
         if ((dir.cx() == 1 && dir.cy() == 0 && dir.cz() == -1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == -1) ||
             (dir.cx() == 1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // top south
         if ((dir.cx() == 0 && dir.cy() == -1 && dir.cz() == 1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == 1) ||
             (dir.cx() == 0 && dir.cy() == -1 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // top north
         if ((dir.cx() == 0 && dir.cy() == 1 && dir.cz() == 1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == 1) ||
             (dir.cx() == 0 && dir.cy() == 1 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(1)),
                                                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // top west
         if ((dir.cx() == -1 && dir.cy() == 0 && dir.cz() == 1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == 1) ||
             (dir.cx() == -1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // top east
         if ((dir.cx() == 1 && dir.cy() == 0 && dir.cz() == 1) || (dir.cx() == 0 && dir.cy() == 0 && dir.cz() == 1) ||
             (dir.cx() == 1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(1)),
                                                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // south-west
         if ((dir.cx() == -1 && dir.cy() == -1 && dir.cz() == 0) ||
             (dir.cx() == 0 && dir.cy() == -1 && dir.cz() == 0) || (dir.cx() == -1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // south-east
         if ((dir.cx() == 1 && dir.cy() == -1 && dir.cz() == 0) || (dir.cx() == 0 && dir.cy() == -1 && dir.cz() == 0) ||
             (dir.cx() == 1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // north-west
         if ((dir.cx() == -1 && dir.cy() == 1 && dir.cz() == 0) || (dir.cx() == 0 && dir.cy() == 1 && dir.cz() == 0) ||
             (dir.cx() == -1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }

         // north-east
         if ((dir.cx() == 1 && dir.cy() == 1 && dir.cz() == 0) || (dir.cx() == 0 && dir.cy() == 1 && dir.cz() == 0) ||
             (dir.cx() == 1 && dir.cy() == 0 && dir.cz() == 0))
         {
            const real_t intersection = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(1), real_c(0)),
                                                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                                                *normalFieldIt, interfacePointLocation);
            if (intersection > real_c(0))
            {
               boundaryHandling->setFlag(flagInfo.keepInterfaceForWettingFlag, flagFieldIt.x() + dir.cx(),
                                         flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
               continue;
            }
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

} // namespace free_surface
} // namespace walberla