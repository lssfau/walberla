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
//! \file ThermalExpansion.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

class ThermalExpansion
{
public:
   ThermalExpansion(const uint_t numParticleTypes);
   ThermalExpansion(const ThermalExpansion& other) = default;
   ThermalExpansion(ThermalExpansion&& other) = default;
   ThermalExpansion& operator=(const ThermalExpansion& other) = default;
   ThermalExpansion& operator=(ThermalExpansion&& other) = default;

   void setLinearExpansionCoefficient(const size_t type, const real_t& val);
   real_t getLinearExpansionCoefficient(const size_t type) const;

   template <typename Accessor>
   void operator()(const size_t p_idx, Accessor& ac) const;
private:
       uint_t numParticleTypes_;

   std::vector<real_t> linearExpansionCoefficient_ {};
};

ThermalExpansion::ThermalExpansion(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;

   linearExpansionCoefficient_.resize(numParticleTypes, real_t(0));
}


inline void ThermalExpansion::setLinearExpansionCoefficient(const size_t type, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type, numParticleTypes_ );
   linearExpansionCoefficient_[type] = val;
}


inline real_t ThermalExpansion::getLinearExpansionCoefficient(const size_t type) const
{
   WALBERLA_ASSERT_LESS( type, numParticleTypes_ );
   return linearExpansionCoefficient_[type];
}

template <typename Accessor>
inline void ThermalExpansion::operator()(const size_t p_idx,
                                         Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   const auto Tc = ac.getTemperature(p_idx)-real_t(273);
   ac.setRadiusAtTemperature(p_idx, ac.getRadius273K(p_idx) * (real_t(1) + Tc * getLinearExpansionCoefficient(ac.getType(p_idx))));
   ac.setInteractionRadius(p_idx, ac.getRadiusAtTemperature(p_idx));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
