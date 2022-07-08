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
//! \file SmoothingSweep.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Smooth fill levels (used for finite difference curvature computation).
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/BlockDataID.h"

namespace walberla
{
namespace free_surface
{
// forward declaration
template< typename Stencil_T >
class KernelK8;

/***********************************************************************************************************************
 * Smooth fill levels such that interface-neighboring cells get assigned a new fill level. This is required for
 * computing the interface curvature using the finite difference method.
 *
 * The same smoothing kernel is used as in the dissertation of S. Bogner, 2017, i.e., the K8 kernel with support
 * radius 2 from
 * Williams, Kothe and Puckett, "Accuracy and Convergence of Continuum Surface Tension Models", 1998.
 **********************************************************************************************************************/
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class SmoothingSweep
{
 protected:
   using vector_t = typename std::remove_const< typename VectorField_T::value_type >::type;
   using flag_t   = typename std::remove_const< typename FlagField_T::value_type >::type;

 public:
   SmoothingSweep(const BlockDataID& smoothFillFieldID, const ConstBlockDataID& fillFieldID,
                  const ConstBlockDataID& flagFieldID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                  const Set< FlagUID >& obstacleFlagIDSet, bool includeObstacleNeighbors)
      : smoothFillFieldID_(smoothFillFieldID), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFlagIDSet_(obstacleFlagIDSet),
        includeObstacleNeighbors_(includeObstacleNeighbors), smoothingKernel_(KernelK8< Stencil_T >(real_c(2.0)))
   {}

   void operator()(IBlock* const block);

 private:
   BlockDataID smoothFillFieldID_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;

   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   Set< FlagUID > obstacleFlagIDSet_;

   bool includeObstacleNeighbors_;

   KernelK8< Stencil_T > smoothingKernel_;
}; // class SmoothingSweep

/***********************************************************************************************************************
 * K8 kernel from Williams, Kothe and Puckett, "Accuracy and Convergence of Continuum Surface Tension Models", 1998.
 **********************************************************************************************************************/
template< typename Stencil_T >
class KernelK8
{
 public:
   KernelK8(real_t epsilon) : epsilon_(epsilon) { stencilSize_ = uint_c(std::ceil(epsilon_) - real_c(1)); }

   // equation (11) in Williams et al. (normalization constant A=1 here, result must be normalized outside the kernel)
   inline real_t kernelFunction(const Vector3< real_t >& dirVec) const
   {
      const real_t r_sqr   = dirVec.sqrLength();
      const real_t eps_sqr = epsilon_ * epsilon_;

      if (r_sqr < eps_sqr)
      {
         real_t result = real_c(1.0) - (r_sqr / (eps_sqr));
         result        = result * result * result * result;

         return result;
      }
      else { return real_c(0.0); }
   }

   inline real_t getSupportRadius() const { return epsilon_; }
   inline uint_t getStencilSize() const { return static_cast< uint_t >(stencilSize_); }

 private:
   real_t epsilon_;     // support radius of the kernel
   uint_t stencilSize_; // size of the stencil which defines included neighbors in smoothing

}; // class KernelK8

} // namespace free_surface
} // namespace walberla

#include "SmoothingSweep.impl.h"