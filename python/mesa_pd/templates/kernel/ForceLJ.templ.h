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
   std::vector<real_t> {{param}} {};
   {%- endfor %}
};

ForceLJ::ForceLJ(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   {% for param in parameters %}
   {{param}}.resize(numParticleTypes * numParticleTypes, real_t(0));
   {%- endfor %}
}

{% for param in parameters %}
inline void ForceLJ::set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   {{param}}[numParticleTypes_*type1 + type2] = val;
   {{param}}[numParticleTypes_*type2 + type1] = val;
}
{%- endfor %}

{% for param in parameters %}
inline real_t ForceLJ::get{{param | capFirst}}(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( {{param}}[numParticleTypes_*type1 + type2],
                                {{param}}[numParticleTypes_*type2 + type1],
                                "parameter matrix for {{param}} not symmetric!");
   return {{param}}[numParticleTypes_*type1 + type2];
}
{%- endfor %}

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

      // Add normal force at contact point
      addForceAtomic( p_idx, ac, f );
      addForceAtomic( np_idx, ac, -f );
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
