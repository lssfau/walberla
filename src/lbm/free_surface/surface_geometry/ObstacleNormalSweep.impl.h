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
//! \file ObstacleNormalSweep.impl.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute a mean obstacle normal in interface cells near solid boundary cells.
//
//======================================================================================================================

#include "core/math/Vector3.h"

#include "field/iterators/IteratorMacros.h"

#include "stencil/D3Q27.h"

#include "ObstacleNormalSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
void ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T >::operator()(IBlock* const block)
{
   // do nothing if obstacle normal must not be computed anywhere
   if (!computeInInterfaceCells_ && !computeInObstacleCells_) { return; }

   // fetch fields
   VectorField_T* const obstacleNormalField = block->getData< VectorField_T >(obstacleNormalFieldID_);
   const FlagField_T* const flagField       = block->getData< const FlagField_T >(flagFieldID_);

   // two ghost layers are required in the flag field
   WALBERLA_ASSERT_EQUAL(flagField->nrOfGhostLayers(), uint_c(2));

   // get flags
   const flag_t interfaceFlag              = flagField->getFlag(interfaceFlagID_);
   const flag_t liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);
   const flag_t obstacleFlagMask           = flagField->getMask(obstacleFlagIDSet_);

   // include ghost layer because solid cells might be located in the (outermost global) ghost layer
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(obstacleNormalField, uint_c(1), {
      if (!computeInGhostLayer_ && (!flagField->isInInnerPart(Cell(x, y, z)))) { continue; }

      const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

      const bool computeInInterfaceCell = computeInInterfaceCells_ && isPartOfMaskSet(flagFieldPtr, interfaceFlag) &&
                                          isFlagInNeighborhood< Stencil_T >(flagFieldPtr, obstacleFlagMask);

      const bool computeInObstacleCell = computeInObstacleCells_ && isPartOfMaskSet(flagFieldPtr, obstacleFlagMask) &&
                                         isFlagInNeighborhood< Stencil_T >(flagFieldPtr, liquidInterfaceGasFlagMask);

      // IMPORTANT REMARK: do not restrict this algorithm to obstacle cells that are direct neighbors of interface
      // cells; the succeeding ObstacleFillLevelSweep and SmoothingSweep use the values computed here and the latter
      // must work on an at least two-cell neighborhood of interface cells

      vector_t& obstacleNormal = obstacleNormalField->get(x, y, z);

      if (computeInInterfaceCell)
      {
         WALBERLA_ASSERT(!computeInObstacleCell);

         // compute mean obstacle, i.e., mean wall normal in interface cell (see dissertation of S. Donath, 2011,
         // section 6.3.5.2)
         computeObstacleNormalInInterfaceCell(obstacleNormal, flagFieldPtr, obstacleFlagMask);
      }
      else
      {
         if (computeInObstacleCell)
         {
            WALBERLA_ASSERT(!computeInInterfaceCell);

            // compute mean obstacle normal in obstacle cell
            computeObstacleNormalInObstacleCell(obstacleNormal, flagFieldPtr, liquidInterfaceGasFlagMask);
         }
         else
         {
            // set obstacle normal of all other cells to zero
            obstacleNormal.set(real_c(0), real_c(0), real_c(0));
         }
      }

      // normalize mean obstacle normal
      const real_t sqrObstNormal = obstacleNormal.sqrLength();
      if (sqrObstNormal > real_c(0))
      {
         const real_t invlength = -real_c(1) / real_c(std::sqrt(sqrObstNormal));
         obstacleNormal *= invlength;
      }
   }); // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
template< typename FlagFieldIt_T >
void ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T >::computeObstacleNormalInInterfaceCell(
   vector_t& obstacleNormal, const FlagFieldIt_T& flagFieldIt, const flag_t& obstacleFlagMask)
{
   uint_t obstCount = uint_c(0);
   obstacleNormal   = vector_t(real_c(0));

   for (auto i = Stencil_T::beginNoCenter(); i != Stencil_T::end(); ++i)
   {
      // only consider directions in which there is an obstacle cell
      if (isPartOfMaskSet(flagFieldIt.neighbor(*i), obstacleFlagMask))
      {
         obstacleNormal += vector_t(real_c(-i.cx()), real_c(-i.cy()), real_c(-i.cz()));
         ++obstCount;
      }
   }
   obstacleNormal = obstCount > uint_c(0) ? obstacleNormal / real_c(obstCount) : vector_t(real_c(0));
}

template< typename Stencil_T, typename FlagField_T, typename VectorField_T >
template< typename FlagFieldIt_T >
void ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T >::computeObstacleNormalInObstacleCell(
   vector_t& obstacleNormal, const FlagFieldIt_T& flagFieldIt, const flag_t& liquidInterfaceGasFlagMask)
{
   uint_t obstCount = uint_c(0);
   obstacleNormal   = vector_t(real_c(0));

   for (auto i = Stencil_T::beginNoCenter(); i != Stencil_T::end(); ++i)
   {
      // only consider directions in which there is a liquid, interface, or gas cell
      if (isPartOfMaskSet(flagFieldIt.neighbor(*i), liquidInterfaceGasFlagMask))
      {
         obstacleNormal += vector_t(real_c(i.cx()), real_c(i.cy()), real_c(i.cz()));

         ++obstCount;
      }
   }
   obstacleNormal = obstCount > uint_c(0) ? obstacleNormal / real_c(obstCount) : vector_t(real_c(0));
}

} // namespace free_surface
} // namespace walberla
