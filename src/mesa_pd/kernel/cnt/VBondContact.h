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

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
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
 * VBond interaction model
 *
 * Implementation according to:
 *
 * Ostanin, I., Dumitrică, T., Eibl, S., and Rüde, U. (October 4, 2019).
 * "Size-Independent Mechanical Response of Ultrathin Carbon Nanotube Films in Mesoscopic Distinct Element Method Simulations."
 * ASME. J. Appl. Mech. December 2019; 86(12): 121006.
 * https://doi.org/10.1115/1.4044413
 */
class VBondContact
{
public:
   VBondContact( const Vector3<bool>& isPeriodic = Vector3<bool>{false, false, false},
                 const Vector3<int64_t>& maxSegments = Vector3<real_t>{0, 0, 0} )
   : isPeriodic_(isPeriodic)
   , maxSegments_(maxSegments)
   {}

   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor &ac);

   static constexpr real_t S = 142.7_r; // Area of the bond,
   static constexpr real_t E = 1029_r * 0.006242_r; // Conversion from GPa to eV/A^3
   static constexpr real_t G = 459_r * 0.006242_r;
   static constexpr real_t a = 2_r * 6.78_r; // (10,10) CNTS
   static constexpr real_t J = 3480_r;
   static constexpr real_t Jp = 2_r * J;

   // Stiffnesses, equilibrium length etc
   static constexpr real_t B1 = E * S / a;
   static constexpr real_t B2 = 12_r * E * J / a;
   static constexpr real_t B3 = -2_r * E * J / a - G * Jp / (2_r * a);
   static constexpr real_t B4 = G * Jp / a;



   /// Get tensile energy from last contact.
   auto getLastTensileEnergy() const {return tensileEnergy;}
   /// Get shear energy from last contact.
   auto getLastShearEnergy() const {return shearEnergy;}
   /// Get bending energy from last contact.
   auto getLastBendingEnergy() const {return bendingEnergy;}
   /// Get twisting energy from last contact.
   auto getLastTwistingEnergy() const {return twistingEnergy;}
private:
   real_t tensileEnergy;
   real_t shearEnergy;
   real_t bendingEnergy;
   real_t twistingEnergy;

   Vector3<bool> isPeriodic_;
   Vector3<int64_t> maxSegments_;
};


template<typename Accessor>
inline
void VBondContact::operator()(const size_t p_idx1,
                              const size_t p_idx2,
                              Accessor &ac)
{
   Vec3 ri = ac.getPosition(p_idx1);
   Vec3 rj = ac.getPosition(p_idx2);

   // Fix for the issue of segment's undefined sides
   real_t sign = ac.getSegmentID(p_idx1) <= ac.getSegmentID(p_idx2) ? 1_r : -1_r;
   if (((isPeriodic_[0]) && (std::abs(ac.getSegmentID(p_idx2) - ac.getSegmentID(p_idx1)) == maxSegments_[0] - 1)) ||
       ((isPeriodic_[1]) && (std::abs(ac.getSegmentID(p_idx2) - ac.getSegmentID(p_idx1)) == maxSegments_[1] - 1)))
      sign = -sign; // special case for periodic fibers

   Vec3 ni1 = ac.getRotation(p_idx1).getMatrix() * Vec3(sign, 0_r, 0_r);
   Vec3 ni2 = ac.getRotation(p_idx1).getMatrix() * Vec3(0_r, 1_r, 0_r);
   Vec3 ni3 = ac.getRotation(p_idx1).getMatrix() * Vec3(0_r, 0_r, 1_r);

   Vec3 nj1 = ac.getRotation(p_idx2).getMatrix() * Vec3(-sign, 0_r, 0_r);
   Vec3 nj2 = ac.getRotation(p_idx2).getMatrix() * Vec3(0_r, 1_r, 0_r);
   Vec3 nj3 = ac.getRotation(p_idx2).getMatrix() * Vec3(0_r, 0_r, 1_r);

   // Vectors Dij and dij
   Vec3 Dij = rj - ri;
   real_t s = Dij.length();
   Vec3 dij = Dij * (1_r / s);

   real_t C = 1_r;

   tensileEnergy  = 0.5_r * B1 * (s - a) * (s - a);
   shearEnergy    = B2 * (0.5_r * (nj1 - ni1) * dij - 0.25_r * ni1 * nj1 + 0.75_r);
   bendingEnergy  = (0.25_r * B2 + B3 + 0.5_r * B4) * (ni1 * nj1 + 1_r);
   twistingEnergy = -0.5_r * B4 * (ni1 * nj1 + ni2 * nj2 + ni3 * nj3 - 1_r);

   Vec3 rij = dij;

   Vec3 Fij = C * B1 * (s - a) * rij + B2 / (2_r * s) * ((nj1 - ni1) - ((nj1 - ni1) * dij) * dij);
   Vec3 Fji = -Fij;
   Vec3 M_TB = C * B3 * (nj1 % ni1) - 0.5_r * B4 * (nj2 % ni2 + nj3 % ni3);
   Vec3 Mij = -C * 0.5_r * B2 * (dij % ni1) + M_TB;
   Vec3 Mji =  C * 0.5_r * B2 * (dij % nj1) - M_TB;

   addForceAtomic(p_idx1, ac, Fij);
   addForceAtomic(p_idx2, ac, Fji);
   addTorqueAtomic(p_idx1, ac, Mij);
   addTorqueAtomic(p_idx2, ac, Mji);
}

} //namespace cnt
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla