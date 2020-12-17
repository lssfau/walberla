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

#include <mesa_pd/data/Flags.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

/**
 * Clear all ghost particles and reset ghost owner information
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
 * \post All ghost particles are deleted.
 * \post All ghost owners are reset.
 * \post All cached information is reset.
 *
 * \ingroup mesa_pd_mpi
 */
class ClearGhostOwnerSync
{
public:
   template <typename Accessor>
   void operator()(Accessor& ac) const;
};

template <typename Accessor>
void ClearGhostOwnerSync::operator()(Accessor& ac) const
{
   for (size_t idx = 0; idx < ac.size(); )
   {
      if (data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::GHOST))
      {
         //ghost particle
         idx = ac.erase(idx);
         continue;
      } else
      {
         //local particle
         ac.getGhostOwnersRef(idx).clear();
         ac.getNeighborStateRef(idx).clear();
      }
      ++idx;
   }
}

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla
