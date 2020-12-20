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
//! \author Grigorii Drozdov, <drozd013@umn.edu>
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
 * PFC style damping
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
 *
 * \ingroup mesa_pd_kernel
 */
class PFCDamping
{
public:
   PFCDamping(const real_t alpha) : alpha_(alpha) {}

   template <typename Accessor>
   void operator()(const size_t p_idx, Accessor& ac) const;
private:
   real_t alpha_ = 0_r;
};

template <typename Accessor>
inline void PFCDamping::operator()(const size_t p_idx,
                                   Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   Vec3 damp_F(0,0,0);
   Vec3 damp_M(0,0,0);

   for (size_t i = 0; i < 3; i++)
   {
      damp_F[i] = - alpha_ * std::fabs( ac.getForce(p_idx)[i] ) * math::sign( ac.getLinearVelocity(p_idx)[i] );
      damp_M[i] = - alpha_ * std::fabs( ac.getTorque(p_idx)[i] ) * math::sign( ac.getAngularVelocity(p_idx)[i] );
   }

   ac.setForce (p_idx, ac.getForce(p_idx)  + damp_F);
   ac.setTorque(p_idx, ac.getTorque(p_idx) + damp_M);
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
