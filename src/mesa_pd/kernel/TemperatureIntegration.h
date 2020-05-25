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
//! \file TemperatureIntegration.h
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

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Kernel which explicitly integrates a particle in time.
 * The heat flux is converted into a temperature change.
 *
 * This kernel requires the following particle accessor interface
 * \code
 * const walberla::real_t& getTemperature(const size_t p_idx) const;
 * void setTemperature(const size_t p_idx, const walberla::real_t& v);
 *
 * const walberla::real_t& getHeatFlux(const size_t p_idx) const;
 * void setHeatFlux(const size_t p_idx, const walberla::real_t& v);
 *
 * const walberla::real_t& getInvMass(const size_t p_idx) const;
 *
 * const uint_t& getType(const size_t p_idx) const;
 *
 * \endcode
 *
 * \pre  Heat flux has to be set/reduced.
 * \post Heat flux is reset to 0.
 * \ingroup mesa_pd_kernel
 */
class TemperatureIntegration
{
public:
   TemperatureIntegration(const real_t dt, const uint_t numParticleTypes);
   TemperatureIntegration(const TemperatureIntegration& other) = default;
   TemperatureIntegration(TemperatureIntegration&& other) = default;
   TemperatureIntegration& operator=(const TemperatureIntegration& other) = default;
   TemperatureIntegration& operator=(TemperatureIntegration&& other) = default;

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;

   
   /// assumes this parameter is symmetric
   void setInvSpecificHeat(const size_t type, const real_t& val);

   
   real_t getInvSpecificHeat(const size_t type) const;
private:
   real_t dt_ = real_t(0.0);

   uint_t numParticleTypes_;
   
   std::vector<real_t> invSpecificHeat_ {};
};

TemperatureIntegration::TemperatureIntegration(const real_t dt, const uint_t numParticleTypes)
   : dt_(dt)
{
   numParticleTypes_ = numParticleTypes;
   
   invSpecificHeat_.resize(numParticleTypes, real_t(0));
}


inline void TemperatureIntegration::setInvSpecificHeat(const size_t type, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type, numParticleTypes_ );
   invSpecificHeat_[type] = val;
}


inline real_t TemperatureIntegration::getInvSpecificHeat(const size_t type) const
{
   WALBERLA_ASSERT_LESS( type, numParticleTypes_ );
   return invSpecificHeat_[type];
}

template <typename Accessor>
inline void TemperatureIntegration::operator()(const size_t p_idx,
                                               Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   //formula for heat capacity
   ac.setTemperature(p_idx, getInvSpecificHeat(ac.getType(p_idx)) * ac.getInvMass(p_idx) * ac.getHeatFlux(p_idx) * dt_ + ac.getTemperature(p_idx));
   ac.setHeatFlux   (p_idx, real_t(0));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla