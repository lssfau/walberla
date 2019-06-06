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
//! \file ForceLJ.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Kernel which calculates the Lennard Jones froce between two particles.
 *
 * This kernel uses the type property of a particle to decide on the material parameters.
 *
 * This kernel requires the following particle accessor interface
 * \code
 * const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const;
 *
 * walberla::mesa_pd::Vec3& getForceRef(const size_t p_idx);
 *
 * const uint_t& getType(const size_t p_idx) const;
 *
 * \endcode
 * \ingroup mesa_pd_kernel
 */
class ForceLJ
{
public:
   ForceLJ(const uint_t numParticleTypes);
   ForceLJ(const ForceLJ& other) = default;
   ForceLJ(ForceLJ&& other) = default;
   ForceLJ& operator=(const ForceLJ& other) = default;
   ForceLJ& operator=(ForceLJ&& other) = default;

   template <typename Accessor>
   void operator()(const size_t p_idx, const size_t np_idx, Accessor& ac) const;

   
   /// assumes this parameter is symmetric
   void setEpsilon(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setSigma(const size_t type1, const size_t type2, const real_t& val);

   
   real_t getEpsilon(const size_t type1, const size_t type2) const;
   real_t getSigma(const size_t type1, const size_t type2) const;
private:
   uint_t numParticleTypes_;
   
   std::vector<real_t> epsilon {};
   std::vector<real_t> sigma {};
};

ForceLJ::ForceLJ(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   
   epsilon.resize(numParticleTypes * numParticleTypes, real_t(0));
   sigma.resize(numParticleTypes * numParticleTypes, real_t(0));
}


inline void ForceLJ::setEpsilon(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   epsilon[numParticleTypes_*type1 + type2] = val;
   epsilon[numParticleTypes_*type2 + type1] = val;
}
inline void ForceLJ::setSigma(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   sigma[numParticleTypes_*type1 + type2] = val;
   sigma[numParticleTypes_*type2 + type1] = val;
}


inline real_t ForceLJ::getEpsilon(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( epsilon[numParticleTypes_*type1 + type2],
                                epsilon[numParticleTypes_*type2 + type1],
                                "parameter matrix for epsilon not symmetric!");
   return epsilon[numParticleTypes_*type1 + type2];
}
inline real_t ForceLJ::getSigma(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( sigma[numParticleTypes_*type1 + type2],
                                sigma[numParticleTypes_*type2 + type1],
                                "parameter matrix for sigma not symmetric!");
   return sigma[numParticleTypes_*type1 + type2];
}

template <typename Accessor>
inline void ForceLJ::operator()(const size_t p_idx, const size_t np_idx, Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (p_idx != np_idx)
   {
      Vec3 dir = ac.getPosition(p_idx) - ac.getPosition(np_idx);
      const real_t rsq = sqrLength(dir);
      const real_t sr2 = real_t(1.0) / rsq;
      const real_t sr2sigma = sr2 * getSigma(ac.getType(p_idx), ac.getType(np_idx)) * getSigma(ac.getType(p_idx), ac.getType(np_idx));
      const real_t sr6 = sr2sigma * sr2sigma * sr2sigma;
      const real_t force = real_t(48) * sr6 * ( sr6 - real_t(0.5) ) * sr2 * getEpsilon(ac.getType(p_idx), ac.getType(np_idx));
      const Vec3 f = force * dir;

#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(p_idx)[0]  += f[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(p_idx)[1]  += f[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(p_idx)[2]  += f[2];

#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(np_idx)[0]  -= f[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(np_idx)[1]  -= f[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getForceRef(np_idx)[2]  -= f[2];
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla