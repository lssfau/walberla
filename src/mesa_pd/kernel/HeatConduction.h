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
//! \file HeatConduction.h
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
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Basic DEM kernel
 *
 * This DEM kernel supports spring&dashpot in normal direction as well as friction in tangential direction.
 *
 * \code
 * const walberla::real_t& getTemperature(const size_t p_idx) const;
 *
 * const walberla::real_t& getHeatFlux(const size_t p_idx) const;
 * void setHeatFlux(const size_t p_idx, const walberla::real_t& v);
 * walberla::real_t& getHeatFluxRef(const size_t p_idx);
 *
 * const uint_t& getType(const size_t p_idx) const;
 *
 * \endcode
 * \ingroup mesa_pd_kernel
 */
class HeatConduction
{
public:
   HeatConduction(const uint_t numParticleTypes);
   HeatConduction(const HeatConduction& other) = default;
   HeatConduction(HeatConduction&& other) = default;
   HeatConduction& operator=(const HeatConduction& other) = default;
   HeatConduction& operator=(HeatConduction&& other) = default;

   template <typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor& ac) const;

   
   /// assumes this parameter is symmetric
   void setConductance(const size_t type1, const size_t type2, const real_t& val);

   
   real_t getConductance(const size_t type1, const size_t type2) const;

private:
   uint_t numParticleTypes_;
   
   std::vector<real_t> conductance_ {};
};

HeatConduction::HeatConduction(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   
   conductance_.resize(numParticleTypes * numParticleTypes, real_t(0));
}


inline void HeatConduction::setConductance(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   conductance_[numParticleTypes_*type1 + type2] = val;
   conductance_[numParticleTypes_*type2 + type1] = val;
}


inline real_t HeatConduction::getConductance(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( conductance_[numParticleTypes_*type1 + type2],
                                conductance_[numParticleTypes_*type2 + type1],
                                "parameter matrix for conductance not symmetric!");
   return conductance_[numParticleTypes_*type1 + type2];
}

template <typename Accessor>
inline void HeatConduction::operator()(const size_t p_idx1,
                                       const size_t p_idx2,
                                       Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (p_idx1 != p_idx2)
   {
      auto deltaT = ac.getTemperature(p_idx1) - ac.getTemperature(p_idx2);
      ac.getHeatFluxRef(p_idx1) -= deltaT * getConductance(ac.getType(p_idx1), ac.getType(p_idx2));
      ac.getHeatFluxRef(p_idx2) += deltaT * getConductance(ac.getType(p_idx1), ac.getType(p_idx2));
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla