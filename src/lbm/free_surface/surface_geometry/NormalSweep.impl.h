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
//! \file NormalSweep.impl.h
//! \ingroup surface_geometry
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute interface normal.
//
//======================================================================================================================

#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "field/iterators/IteratorMacros.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include <type_traits>

#include "NormalSweep.h"
namespace walberla
{
namespace free_surface
{

template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(IBlock* const block)
{
   // fetch fields
   VectorField_T* const normalField     = block->getData< VectorField_T >(normalFieldID_);
   const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);
   const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);

   // two ghost layers are required in the flag field
   WALBERLA_ASSERT_EQUAL(flagField->nrOfGhostLayers(), uint_c(2));

   // get flags
   const flag_t interfaceFlag              = flagField->getFlag(interfaceFlagID_);
   const flag_t liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);
   const flag_t obstacleFlagMask           = flagField->getMask(obstacleFlagIDSet_);

   // evaluate flags in D2Q19 neighborhood for 2D, and in D3Q27 neighborhood for 3D simulations
   using NeighborhoodStencil_T =
      typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;

   // include ghost layer because solid cells might be located in the (outermost global) ghost layer
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(normalField, uint_c(1), {
      if (!computeInGhostLayer_ && (!flagField->isInInnerPart(Cell(x, y, z)))) { continue; }

      const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);
      const typename ScalarField_T::ConstPtr fillFieldPtr(*fillField, x, y, z);

      const bool computeNormalInCell =
         isFlagSet(flagFieldPtr, interfaceFlag) ||
         (computeInInterfaceNeighbors_ && isFlagInNeighborhood< NeighborhoodStencil_T >(flagFieldPtr, interfaceFlag));

      vector_t& normal = normalField->get(x, y, z);

      if (computeNormalInCell)
      {
         if (includeObstacleNeighbors_)
         {
            // requires meaningful fill level values in obstacle cells, as set by ObstacleFillLevelSweep when using
            // curvature computation via the finite difference method
            normal_computation::computeNormal< Stencil_T >(normal, fillFieldPtr, flagFieldPtr,
                                                           liquidInterfaceGasFlagMask | obstacleFlagMask);
         }
         else
         {
            if (modifyNearObstacles_ && isFlagInNeighborhood< Stencil_T >(flagFieldPtr, obstacleFlagMask))
            {
               // near solid boundary cells, compute a Parker-Youngs approximation around a virtual (constructed)
               // midpoint that is displaced by a distance of 0.5 away from the boundary (see dissertation of S. Donath,
               // 2011, section 6.3.5.1); use only for curvature computation based on local triangulation
               normal_computation::computeNormalNearSolidBoundary< Stencil_T >(
                  normal, fillFieldPtr, flagFieldPtr, liquidInterfaceGasFlagMask, obstacleFlagMask);
            }
            else
            {
               normal_computation::computeNormal< Stencil_T >(normal, fillFieldPtr, flagFieldPtr,
                                                              liquidInterfaceGasFlagMask);
            }
         }

         // normalize and negate normal (to make it point from liquid to gas)
         normal = real_c(-1) * normal.getNormalizedOrZero();
      }
      else { normal.set(real_c(0), real_c(0), real_c(0)); }
   }); // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

namespace normal_computation
{
template< typename Stencil_T, typename vector_t, typename ScalarFieldIt_T, typename FlagFieldIt_T, typename flag_t >
void computeNormal(vector_t& normal, const ScalarFieldIt_T& fillFieldIt, const FlagFieldIt_T& flagFieldIt,
                   const flag_t& validNeighborFlagMask)
{
   // All computations are performed in double precision here (and truncated later). This is done to avoid an issue
   // observed with the Intel 19 compiler when built in "DebugOptimized" mode with single precision. There, the result
   // is dependent on the order of the "tmp_*" variables' definitions. It is assumed that an Intel-specific optimization
   // leads to floating point inaccuracies.
   normal = vector_t(real_c(0));

   // avoid accessing neighbors that are out-of-range, i.e., restrict neighbor access to first ghost layer
   const bool useW = flagFieldIt.x() >= cell_idx_c(0);
   const bool useE = flagFieldIt.x() < cell_idx_c(flagFieldIt.getField()->xSize());
   const bool useS = flagFieldIt.y() >= cell_idx_c(0);
   const bool useN = flagFieldIt.y() < cell_idx_c(flagFieldIt.getField()->ySize());

   // loops are unrolled for improved computational performance
   // IMPORTANT REMARK: the non-unrolled implementation was observed to give different results at O(1e-15); this
   // accumulated and lead to inaccuracies, e.g., a drop wetting a surface became asymmetrical and started to move
   // sideways
   if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 >)
   {
      // get fill level in neighboring cells
      const double tmp_S  = useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(0, -1, 0)) :
                               static_cast< double >(0);
      const double tmp_N  = useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(0, 1, 0)) :
                               static_cast< double >(0);
      const double tmp_W  = useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(-1, 0, 0)) :
                               static_cast< double >(0);
      const double tmp_E  = useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(1, 0, 0)) :
                               static_cast< double >(0);
      const double tmp_SW = useS && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, -1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(-1, -1, 0)) :
                               static_cast< double >(0);
      const double tmp_SE = useS && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, -1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(1, -1, 0)) :
                               static_cast< double >(0);
      const double tmp_NW = useN && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(-1, 1, 0)) :
                               static_cast< double >(0);
      const double tmp_NE = useN && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 1, 0), validNeighborFlagMask) ?
                               static_cast< double >(fillFieldIt.neighbor(1, 1, 0)) :
                               static_cast< double >(0);

      // compute normal with Parker-Youngs approximation (PY)
      const double weight_1 = real_c(4);
      normal[0]             = real_c(weight_1 * (tmp_E - tmp_W));
      normal[1]             = real_c(weight_1 * (tmp_N - tmp_S));
      normal[2]             = real_c(0);

      const double weight_2 = real_c(2);
      normal[0] += real_c(weight_2 * ((tmp_NE + tmp_SE) - (tmp_NW + tmp_SW)));
      normal[1] += real_c(weight_2 * ((tmp_NE + tmp_NW) - (tmp_SE + tmp_SW)));
   }
   else
   {
      const bool useB = flagFieldIt.z() >= cell_idx_c(0);
      const bool useT = flagFieldIt.z() < cell_idx_c(flagFieldIt.getField()->zSize());

      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
      {
         // get fill level in neighboring cells
         const double tmp_S  = useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, -1, 0)) :
                                  static_cast< double >(0);
         const double tmp_N  = useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, 1, 0)) :
                                  static_cast< double >(0);
         const double tmp_W  = useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(-1, 0, 0)) :
                                  static_cast< double >(0);
         const double tmp_E  = useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(1, 0, 0)) :
                                  static_cast< double >(0);
         const double tmp_SW = useS && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, -1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(-1, -1, 0)) :
                                  static_cast< double >(0);
         const double tmp_SE = useS && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, -1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(1, -1, 0)) :
                                  static_cast< double >(0);
         const double tmp_NW = useN && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(-1, 1, 0)) :
                                  static_cast< double >(0);
         const double tmp_NE = useN && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 1, 0), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(1, 1, 0)) :
                                  static_cast< double >(0);
         const double tmp_B  = useB && isPartOfMaskSet(flagFieldIt.neighbor(0, 0, -1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, 0, -1)) :
                                  static_cast< double >(0);
         const double tmp_T  = useT && isPartOfMaskSet(flagFieldIt.neighbor(0, 0, 1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, 0, 1)) :
                                  static_cast< double >(0);
         const double tmp_BS = useB && useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, -1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, -1, -1)) :
                                  static_cast< double >(0);
         const double tmp_BN = useB && useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, -1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, 1, -1)) :
                                  static_cast< double >(0);
         const double tmp_BW = useB && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, -1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(-1, 0, -1)) :
                                  static_cast< double >(0);
         const double tmp_BE = useB && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, -1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(1, 0, -1)) :
                                  static_cast< double >(0);
         const double tmp_TS = useT && useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, 1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, -1, 1)) :
                                  static_cast< double >(0);
         const double tmp_TN = useT && useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, 1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(0, 1, 1)) :
                                  static_cast< double >(0);
         const double tmp_TW = useT && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, 1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(-1, 0, 1)) :
                                  static_cast< double >(0);
         const double tmp_TE = useT && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, 1), validNeighborFlagMask) ?
                                  static_cast< double >(fillFieldIt.neighbor(1, 0, 1)) :
                                  static_cast< double >(0);

         // compute normal with Parker-Youngs approximation (PY)
         const double weight_1 = real_c(4);
         normal[0]             = real_c(weight_1 * (tmp_E - tmp_W));
         normal[1]             = real_c(weight_1 * (tmp_N - tmp_S));
         normal[2]             = real_c(weight_1 * (tmp_T - tmp_B));

         const double weight_2 = real_c(2);
         normal[0] += real_c(weight_2 * ((tmp_NE + tmp_SE + tmp_TE + tmp_BE) - (tmp_NW + tmp_SW + tmp_TW + tmp_BW)));
         normal[1] += real_c(weight_2 * ((tmp_NE + tmp_NW + tmp_TN + tmp_BN) - (tmp_SE + tmp_SW + tmp_TS + tmp_BS)));
         normal[2] += real_c(weight_2 * ((tmp_TN + tmp_TS + tmp_TE + tmp_TW) - (tmp_BN + tmp_BS + tmp_BE + tmp_BW)));
      }
      else
      {
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >)
         {
            // get fill level in neighboring cells
            const double tmp_S = useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, 0), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(0, -1, 0)) :
                                    static_cast< double >(0);
            const double tmp_N = useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, 0), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(0, 1, 0)) :
                                    static_cast< double >(0);
            const double tmp_W = useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, 0), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(-1, 0, 0)) :
                                    static_cast< double >(0);
            const double tmp_E = useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, 0), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(1, 0, 0)) :
                                    static_cast< double >(0);
            const double tmp_SW =
               useS && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, -1, 0), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, -1, 0)) :
                  static_cast< double >(0);
            const double tmp_SE =
               useS && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, -1, 0), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, -1, 0)) :
                  static_cast< double >(0);
            const double tmp_NW =
               useN && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 1, 0), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, 1, 0)) :
                  static_cast< double >(0);
            const double tmp_NE =
               useN && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 1, 0), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, 1, 0)) :
                  static_cast< double >(0);
            const double tmp_B = useB && isPartOfMaskSet(flagFieldIt.neighbor(0, 0, -1), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(0, 0, -1)) :
                                    static_cast< double >(0);
            const double tmp_T = useT && isPartOfMaskSet(flagFieldIt.neighbor(0, 0, 1), validNeighborFlagMask) ?
                                    static_cast< double >(fillFieldIt.neighbor(0, 0, 1)) :
                                    static_cast< double >(0);
            const double tmp_BS =
               useB && useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(0, -1, -1)) :
                  static_cast< double >(0);
            const double tmp_BN =
               useB && useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(0, 1, -1)) :
                  static_cast< double >(0);
            const double tmp_BW =
               useB && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, 0, -1)) :
                  static_cast< double >(0);
            const double tmp_BE =
               useB && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, 0, -1)) :
                  static_cast< double >(0);
            const double tmp_TS =
               useT && useS && isPartOfMaskSet(flagFieldIt.neighbor(0, -1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(0, -1, 1)) :
                  static_cast< double >(0);
            const double tmp_TN =
               useT && useN && isPartOfMaskSet(flagFieldIt.neighbor(0, 1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(0, 1, 1)) :
                  static_cast< double >(0);
            const double tmp_TW =
               useT && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 0, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, 0, 1)) :
                  static_cast< double >(0);
            const double tmp_TE =
               useT && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 0, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, 0, 1)) :
                  static_cast< double >(0);
            const double tmp_BSW =
               useB && useS && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, -1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, -1, -1)) :
                  static_cast< double >(0);
            const double tmp_BNW =
               useB && useN && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, 1, -1)) :
                  static_cast< double >(0);
            const double tmp_BSE =
               useB && useS && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, -1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, -1, -1)) :
                  static_cast< double >(0);
            const double tmp_BNE =
               useB && useN && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 1, -1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, 1, -1)) :
                  static_cast< double >(0);
            const double tmp_TSW =
               useT && useS && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, -1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, -1, 1)) :
                  static_cast< double >(0);
            const double tmp_TNW =
               useT && useN && useW && isPartOfMaskSet(flagFieldIt.neighbor(-1, 1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(-1, 1, 1)) :
                  static_cast< double >(0);
            const double tmp_TSE =
               useT && useS && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, -1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, -1, 1)) :
                  static_cast< double >(0);
            const double tmp_TNE =
               useT && useN && useE && isPartOfMaskSet(flagFieldIt.neighbor(1, 1, 1), validNeighborFlagMask) ?
                  static_cast< double >(fillFieldIt.neighbor(1, 1, 1)) :
                  static_cast< double >(0);

            // compute normal with Parker-Youngs approximation (PY)
            const double weight_1 = real_c(4);
            normal[0]             = real_c(weight_1 * (tmp_E - tmp_W));
            normal[1]             = real_c(weight_1 * (tmp_N - tmp_S));
            normal[2]             = real_c(weight_1 * (tmp_T - tmp_B));

            const double weight_2 = real_c(2);
            normal[0] += real_c(weight_2 * ((tmp_NE + tmp_SE + tmp_TE + tmp_BE) - (tmp_NW + tmp_SW + tmp_TW + tmp_BW)));
            normal[1] += real_c(weight_2 * ((tmp_NE + tmp_NW + tmp_TN + tmp_BN) - (tmp_SE + tmp_SW + tmp_TS + tmp_BS)));
            normal[2] += real_c(weight_2 * ((tmp_TN + tmp_TS + tmp_TE + tmp_TW) - (tmp_BN + tmp_BS + tmp_BE + tmp_BW)));

            // weight (=1) corresponding to Parker-Youngs approximation
            normal[0] += real_c((tmp_TNE + tmp_TSE + tmp_BNE + tmp_BSE) - (tmp_TNW + tmp_TSW + tmp_BNW + tmp_BSW));
            normal[1] += real_c((tmp_TNE + tmp_TNW + tmp_BNE + tmp_BNW) - (tmp_TSE + tmp_TSW + tmp_BSE + tmp_BSW));
            normal[2] += real_c((tmp_TNE + tmp_TNW + tmp_TSE + tmp_TSW) - (tmp_BNE + tmp_BNW + tmp_BSE + tmp_BSW));
         }
         else { WALBERLA_ABORT("The chosen stencil type is not implemented in computeNormal()."); }
      }
   }
}

template< typename Stencil_T, typename vector_t, typename ScalarFieldIt_T, typename FlagFieldIt_T, typename flag_t >
void computeNormalNearSolidBoundary(vector_t& normal, const ScalarFieldIt_T& fillFieldIt,
                                    const FlagFieldIt_T& flagFieldIt, const flag_t& validNeighborFlagMask,
                                    const flag_t& obstacleFlagMask)
{
   Vector3< real_t > midPoint(real_c(0));

   // construct the virtual midpoint
   for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
   {
      if (isPartOfMaskSet(flagFieldIt.neighbor(*dir), validNeighborFlagMask) &&
          isPartOfMaskSet(flagFieldIt.neighbor(dir.inverseDir()), obstacleFlagMask))
      {
         if constexpr (Stencil_T::D == uint_t(2))
         {
            midPoint[0] += real_c(dir.cx()) *
                           real_c(stencil::gaussianMultipliers[stencil::D3Q27::idx[stencil::map2Dto3D[2][*dir]]]);
            midPoint[1] += real_c(dir.cy()) *
                           real_c(stencil::gaussianMultipliers[stencil::D3Q27::idx[stencil::map2Dto3D[2][*dir]]]);
            midPoint[2] += real_c(0);
         }
         else
         {
            midPoint[0] += real_c(dir.cx()) * real_c(stencil::gaussianMultipliers[dir.toIdx()]);
            midPoint[1] += real_c(dir.cy()) * real_c(stencil::gaussianMultipliers[dir.toIdx()]);
            midPoint[2] += real_c(dir.cz()) * real_c(stencil::gaussianMultipliers[dir.toIdx()]);
         }
      }
   }

   // restrict the displacement of the virtual midpoint to an absolute value of 0.5
   for (uint_t i = uint_c(0); i != uint_c(3); ++i)
   {
      if (midPoint[i] > real_c(0.0)) { midPoint[i] = real_c(0.5); }
      else
      {
         if (midPoint[i] < real_c(0.0)) { midPoint[i] = real_c(-0.5); }
         // else midPoint[i] == 0
      }
   }

   normal.set(real_c(0), real_c(0), real_c(0));

   // use standard Parker-Youngs approximation (PY) for all cells without solid boundary neighbors
   // otherwise shift neighboring cell by virtual midpoint (also referred to as narrower PY)
   for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
   {
      // skip directions that have obstacleFlagMask set in this and the opposing direction; this is NOT documented in
      // literature, however, it avoids that an undefined fill level is used from the wall cells
      if (isPartOfMaskSet(flagFieldIt.neighbor(*dir), obstacleFlagMask) &&
          isPartOfMaskSet(flagFieldIt.neighbor(dir.inverseDir()), obstacleFlagMask))
      {
         continue;
      }

      cell_idx_t modCx = dir.cx();
      cell_idx_t modCy = dir.cy();
      cell_idx_t modCz = dir.cz();

      // shift neighboring cells by midpoint if they are solid boundary cells
      if (isPartOfMaskSet(flagFieldIt.neighbor(*dir), obstacleFlagMask) ||
          isPartOfMaskSet(flagFieldIt.neighbor(dir.inverseDir()), obstacleFlagMask))
      {
         // truncate cells towards 0, i.e., make the access pattern narrower
         modCx = cell_idx_c(real_c(modCx) + midPoint[0]);
         modCy = cell_idx_c(real_c(modCy) + midPoint[1]);
         modCz = cell_idx_c(real_c(modCz) + midPoint[2]);
      }

      real_t fill;

      if (isPartOfMaskSet(flagFieldIt.neighbor(modCx, modCy, modCz), validNeighborFlagMask))
      {
         // compute normal with formula from regular Parker-Youngs approximation
         if constexpr (Stencil_T::D == uint_t(2))
         {
            fill = fillFieldIt.neighbor(modCx, modCy, modCz) *
                   real_c(stencil::gaussianMultipliers[stencil::D3Q27::idx[stencil::map2Dto3D[2][*dir]]]);
         }
         else { fill = fillFieldIt.neighbor(modCx, modCy, modCz) * real_c(stencil::gaussianMultipliers[dir.toIdx()]); }

         normal[0] += real_c(dir.cx()) * fill;
         normal[1] += real_c(dir.cy()) * fill;
         normal[2] += real_c(dir.cz()) * fill;
      }
      else { normal = Vector3< real_t >(real_c(0)); }
   }
}
} // namespace normal_computation
} // namespace free_surface
} // namespace walberla
