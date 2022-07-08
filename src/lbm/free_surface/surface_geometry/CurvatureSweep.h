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
//! \file CurvatureSweep.h
//! \ingroup surface_geometry
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweeps for computing the interface curvature.
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"
#include "stencil/Directions.h"

#include <type_traits>
#include <vector>

#include "ContactAngle.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Compute the interface curvature using a finite difference scheme (Parker-Youngs approach) as described in
 * - dissertation of S. Bogner, 2017 (section 4.4.2.1)
 * which is based on
 * - Brackbill, Kothe and Zemach, "A continuum method ...", 1992
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class CurvatureSweepFiniteDifferences
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   CurvatureSweepFiniteDifferences(const BlockDataID& curvatureFieldID, const ConstBlockDataID& normalFieldID,
                                   const ConstBlockDataID& obstacleNormalFieldID, const ConstBlockDataID& flagFieldID,
                                   const FlagUID& interfaceFlagID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                                   const Set< FlagUID >& obstacleFlagIDSet, bool enableWetting,
                                   const ContactAngle& contactAngle)
      : curvatureFieldID_(curvatureFieldID), normalFieldID_(normalFieldID),
        obstacleNormalFieldID_(obstacleNormalFieldID), flagFieldID_(flagFieldID), enableWetting_(enableWetting),
        contactAngle_(contactAngle), interfaceFlagID_(interfaceFlagID),
        liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFlagIDSet_(obstacleFlagIDSet)
   {}

   void operator()(IBlock* const block);

   /********************************************************************************************************************
    * Returns an adjusted interface normal according to the wetting model from the dissertation of S. Bogner, 2017
    * (section 4.4.2.1).
    *******************************************************************************************************************/
   template< typename VectorIt_T >
   Vector3< real_t > getNormalWithWetting(VectorIt_T normalFieldIt, VectorIt_T obstacleNormalFieldIt,
                                          const stencil::Direction dir)
   {
      // get reversed interface normal
      const Vector3< real_t > n = -*normalFieldIt;

      // get n_w, i.e., obstacle normal
      const Vector3< real_t > nw = obstacleNormalFieldIt.neighbor(dir);

      // get n_t: vector tangent to the wall and normal to the contact line; obtained by subtracting the wall-normal
      // component from the interface normal n
      Vector3< real_t > nt = n - (nw * n) * nw; // "(nw * n) * nw" is orthogonal projection of nw on n
      nt                   = nt.getNormalizedOrZero();

      // compute interface normal at wall according to wetting model (equation 4.21 in dissertation of S. Bogner,
      // 2017)
      const Vector3< real_t > nWall = nw * contactAngle_.getCos() + nt * contactAngle_.getSin();

      // extrapolate into obstacle cell to obtain boundary value; expression comes from 1D extrapolation formula
      // IMPORTANT REMARK: corner and diagonal directions must use different formula, Bogner did not consider this in
      // his implementation; here, nevertheless Bogner's approach is used
      const Vector3< real_t > nWallExtrapolated = n + real_c(2) * (nWall - n);

      return -nWallExtrapolated;
   }

 private:
   BlockDataID curvatureFieldID_;
   ConstBlockDataID normalFieldID_;
   ConstBlockDataID obstacleNormalFieldID_;
   ConstBlockDataID flagFieldID_;

   bool enableWetting_;
   ContactAngle contactAngle_;

   FlagUID interfaceFlagID_;
   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   Set< FlagUID > obstacleFlagIDSet_;
}; // class CurvatureSweepFiniteDifferences

/***********************************************************************************************************************
 * Compute the interface curvature using local triangulation as described in
 * - dissertation of T. Pohl, 2008 (section 2.5)
 * - dissertation of S. Donath, 2011 (wetting model, section 6.3.3)
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class CurvatureSweepLocalTriangulation
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;

 public:
   CurvatureSweepLocalTriangulation(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                                    const BlockDataID& curvatureFieldID, const ConstBlockDataID& normalFieldID,
                                    const ConstBlockDataID& fillFieldID, const ConstBlockDataID& flagFieldID,
                                    const ConstBlockDataID& obstacleNormalFieldID, const FlagUID& interfaceFlagID,
                                    const Set< FlagUID >& obstacleFlagIDSet, bool enableWetting,
                                    const ContactAngle& contactAngle)
      : blockForest_(blockForest), curvatureFieldID_(curvatureFieldID), normalFieldID_(normalFieldID),
        fillFieldID_(fillFieldID), flagFieldID_(flagFieldID), obstacleNormalFieldID_(obstacleNormalFieldID),
        enableWetting_(enableWetting), contactAngle_(contactAngle), interfaceFlagID_(interfaceFlagID),
        obstacleFlagIDSet_(obstacleFlagIDSet)
   {
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 >)
      {
         WALBERLA_ABORT(
            "Curvature computation with local triangulation using a D2Q9 stencil has not been thoroughly tested.");
      }
   }

   void operator()(IBlock* const block);

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   BlockDataID curvatureFieldID_;
   ConstBlockDataID normalFieldID_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID obstacleNormalFieldID_;

   bool enableWetting_;
   ContactAngle contactAngle_;

   FlagUID interfaceFlagID_;
   Set< FlagUID > obstacleFlagIDSet_;
}; // class CurvatureSweepLocalTriangulation

/***********************************************************************************************************************
 * Compute the interface curvature with a simplistic finite difference method. This approach is not documented in
 * literature and neither thoroughly tested or validated.
 * Use it with caution and preferably for testing purposes only.
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class CurvatureSweepSimpleFiniteDifferences
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;

 public:
   CurvatureSweepSimpleFiniteDifferences(const BlockDataID& curvatureFieldID, const ConstBlockDataID& normalFieldID,
                                         const ConstBlockDataID& flagFieldID, const FlagUID& interfaceFlagID,
                                         const Set< FlagUID >& obstacleFlagIDSet, bool enableWetting,
                                         const ContactAngle& contactAngle)
      : curvatureFieldID_(curvatureFieldID), normalFieldID_(normalFieldID), flagFieldID_(flagFieldID),
        enableWetting_(enableWetting), contactAngle_(contactAngle), interfaceFlagID_(interfaceFlagID),
        obstacleFlagIDSet_(obstacleFlagIDSet)
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
         "You are using curvature computation based on a simplistic finite difference method. This "
         "was implemented for testing purposes only and has not been thoroughly "
         "validated and tested in the current state of the code. Use it with caution.");
   }

   void operator()(IBlock* const block);

 private:
   BlockDataID curvatureFieldID_;
   ConstBlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;

   bool enableWetting_;
   ContactAngle contactAngle_;

   FlagUID interfaceFlagID_;
   Set< FlagUID > obstacleFlagIDSet_;
}; // class CurvatureSweepSimpleFiniteDifferences

} // namespace free_surface
} // namespace walberla

#include "CurvatureSweep.impl.h"