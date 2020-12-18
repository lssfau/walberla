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
 * Isotropic vdW contact model.
 *
 * implementation follows:
 * I. Ostanin, R. Ballarini, D. Potyondy, T. Dumitrica, A distinct element method for large scale simulations of carbon nanotube assemblies, J. Mech. Phys. Solid. 61 (2013) 762-782.
 * https://doi.org/10.1016/j.jmps.2012.10.016
 */
class IsotropicVDWContact
{
public:
   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor &ac);
   static real_t equilibriumDistance();

   static constexpr real_t eps = 0.07124_r;
   static constexpr real_t A = 0.0223_r;
   static constexpr real_t B = 1.31_r;
   static constexpr real_t alpha = 9.5_r;
   static constexpr real_t beta = 4.0_r;
   static constexpr real_t r = 6.78_r;
   static constexpr real_t Rinv = 1.0_r / r;
   static constexpr real_t Dc = 0.4_r;

   auto getLastForce() const {return F_;}
   auto getLastEnergy() const {return U_;}
private:
   real_t F_ = 0_r; ///< resulting force from the last interaction
   real_t U_ = 0_r; ///< resulting energy from the last interaction
};

template<typename Accessor>
inline
void IsotropicVDWContact::operator()(const size_t p_idx1,
                                     const size_t p_idx2,
                                     Accessor &ac)
{
   Vec3 n = ac.getPosition(p_idx2) - ac.getPosition(p_idx1); ///< contact normal
   real_t L = n.length();
   n *= (1_r/L);

   real_t D = L * Rinv - 2_r;
   F_ = 0.01_r; //default value
   U_ = 0_r; //default value

   if (D > Dc)
   {
      const auto pow_alpha = std::pow(D, alpha);
      const auto pow_beta = std::pow(D, beta);
      F_ = 4_r * eps * (-(alpha * A) / (pow_alpha * D) + (beta * B) / (pow_beta * D));
      U_ = 4_r * eps * (A / pow_alpha - B / pow_beta);
   }

   Vec3 force = n * F_;
   addForceAtomic(p_idx1, ac,  force);
   addForceAtomic(p_idx2, ac, -force);
}

real_t IsotropicVDWContact::equilibriumDistance()
{
   return r * ( std::pow( (alpha*A)/(beta*B), 1_r/(alpha-beta)) + 2_r);
}

} //namespace cnt
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla