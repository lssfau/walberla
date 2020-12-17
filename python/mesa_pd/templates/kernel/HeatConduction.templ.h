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
 * Heat conduction interaction kernel
 *
 * This kernel implements a simple heat conduction: \frac{dQ}{dt} = - \alpha (T_2 - T_1)
 *
 * \code
   {%- for prop in interface %}
   {%- if 'g' in prop.access %}
 * const {{prop.type}}& get{{prop.name | capFirst}}(const size_t p_idx) const;
   {%- endif %}
   {%- if 's' in prop.access %}
 * void set{{prop.name | capFirst}}(const size_t p_idx, const {{prop.type}}& v);
   {%- endif %}
   {%- if 'r' in prop.access %}
 * {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t p_idx);
   {%- endif %}
 *
   {%- endfor %}
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

   {% for param in parameters %}
   /// assumes this parameter is symmetric
   void set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val);
   {%- endfor %}

   {% for param in parameters %}
   real_t get{{param | capFirst}}(const size_t type1, const size_t type2) const;
   {%- endfor %}

private:
   uint_t numParticleTypes_;
   {% for param in parameters %}
   std::vector<real_t> {{param}}_ {};
   {%- endfor %}
};

HeatConduction::HeatConduction(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   {% for param in parameters %}
   {{param}}_.resize(numParticleTypes * numParticleTypes, real_t(0));
   {%- endfor %}
}

{% for param in parameters %}
inline void HeatConduction::set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   {{param}}_[numParticleTypes_*type1 + type2] = val;
   {{param}}_[numParticleTypes_*type2 + type1] = val;
}
{%- endfor %}

{% for param in parameters %}
inline real_t HeatConduction::get{{param | capFirst}}(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( {{param}}_[numParticleTypes_*type1 + type2],
                                {{param}}_[numParticleTypes_*type2 + type1],
                                "parameter matrix for {{param}} not symmetric!");
   return {{param}}_[numParticleTypes_*type1 + type2];
}
{%- endfor %}

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
