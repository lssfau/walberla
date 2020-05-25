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
   ThermalExpansion() = default;

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;
};

template <typename Accessor>
inline void ThermalExpansion::operator()(const size_t idx,
                                         Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   ac.setRadius(idx, real_t(0.004)  + (ac.getTemperature(idx)-real_t(273)) * real_t(0.002) * real_t(0.001));
   ac.setInteractionRadius(idx, ac.getRadius(idx));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
