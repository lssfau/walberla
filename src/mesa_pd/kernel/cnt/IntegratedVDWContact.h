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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <core/math/Constants.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
namespace cnt {

/**
 * vdW Contact with integration
 */
class IntegratedVDWContact
{
public:
   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor &ac);

   static constexpr real_t R_CNT = 6.78_r; ///< CNT radius
   static constexpr real_t T = 2_r * R_CNT; ///< Height of a cylindrical segment
   static constexpr real_t eps = 0.07124_r;
   static constexpr real_t A = 0.0223_r;
   static constexpr real_t B = 1.31_r;
   static constexpr real_t alpha = 9.5_r;
   static constexpr real_t beta = 4.0_r;

   // vdW adhesion + linear repulsion potential.
   const real_t r0 = R_CNT * (2_r + std::pow((alpha * A / (beta * B)), 1_r / (alpha - beta)));
   const real_t u0 = 4_r * eps * (A / std::pow(r0 / R_CNT - 2_r, (alpha)) - B / std::pow(r0 / R_CNT - 2_r, (beta)));

   auto getLastEnergy() const {return energy_;}
private:
   real_t energy_; ///< total potential
};


template<typename Accessor>
inline
void IntegratedVDWContact::operator()(const size_t p_idx1,
                                      const size_t p_idx2,
                                      Accessor &ac)
{
   constexpr real_t K_n = 200.0_r; // Good for fabrics modeling

   // particle centers
   Vec3 O1 = ac.getPosition(p_idx1);
   Vec3 O2 = ac.getPosition(p_idx2);

   // axial vectors
   Vec3 b1 = ac.getRotation(p_idx1).getMatrix() * Vec3(1_r, 0_r, 0_r);
   Vec3 b2 = ac.getRotation(p_idx2).getMatrix() * Vec3(1_r, 0_r, 0_r);

   energy_ = 0_r;
   Vec3 force12(0);  // Force 12
   Vec3 force21(0);  // Force 21
   Vec3 moment12(0); // Total torque 12
   Vec3 moment21(0); // Total torque 21

   constexpr int Np = 5; // Number of integration points over each axis
   constexpr real_t Npinv = 1.0_r / real_t(Np);
   for (int i = 0; i < Np; ++i) // integral dl1
   {
      for (int j = 0; j < Np; ++j) // integral dl2
      {
         // Levers
         Vec3 l1 = (-0.5_r * T + (0.5_r + real_t(i)) * T * Npinv) * b1;
         Vec3 l2 = (-0.5_r * T + (0.5_r + real_t(j)) * T * Npinv) * b2;

         /// radius vector between dl1 and dl2
         Vec3 R12 = (O2 + l2) - (O1 + l1);

         real_t r12 = R12.length();
         Vec3 n12 = R12 * (1_r / r12);

         real_t dU = 0_r;
         Vec3 dforce12(0);

         if (r12 < r0)
         {
            // elastic interaction
            dU = u0 + K_n * (r12 - r0) * (r12 - r0) * Npinv * Npinv;
            dforce12 = n12 * K_n * (r12 - r0) * Npinv * Npinv;
         } else
         {
            // vdW interaction
            const real_t normDistance = r12 / R_CNT - 2_r;
            const real_t powAlpha = std::pow(normDistance, alpha);
            const real_t powBeta = std::pow(normDistance, beta);
            dU = 4_r * eps * (A / powAlpha - B / powBeta) * Npinv * Npinv;
            dforce12 = n12 * 4_r * eps / R_CNT * Npinv * Npinv *
                       (-(alpha * A) / (powAlpha * normDistance) + (beta * B) / (powBeta * normDistance));
         }

         Vec3 dmoment12 = l2 % dforce12;
         Vec3 dmoment21 = l1 % (-dforce12);

         energy_ += dU;
         force12 += dforce12;
         force21 -= dforce12;
         moment12 += dmoment12;
         moment21 += dmoment21;
      }
   }

   addForceAtomic(p_idx1, ac, force12);
   addForceAtomic(p_idx2, ac, force21);
   addTorqueAtomic(p_idx1, ac,  -moment21);
   addTorqueAtomic(p_idx2, ac,  -moment12);
}

} //namespace cnt
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla