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
//! \file SingleCast.h
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
#include <mesa_pd/data/shape/BaseShape.h>
{%- for shape in particle.shapes %}
#include <mesa_pd/data/shape/{{shape}}.h>
{%- endfor %}

#include <core/Abort.h>
#include <core/debug/Debug.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
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
class SingleCast
{
public:
   template <typename Accessor, typename func, typename... Args>
   auto operator()( size_t idx, Accessor& ac, func& f, Args&&... args );
};

template <typename Accessor, typename func, typename... Args>
auto SingleCast::operator()( size_t idx, Accessor& ac, func& f, Args&&... args )
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   using namespace mesa_pd::data;
   switch (ac.getShape(idx)->getShapeType())
   {
      {%- for shape in particle.shapes %}
      case {{shape}}::SHAPE_TYPE : return f(idx, *static_cast<{{shape}}*>(ac.getShape(idx)), std::forward<Args>(args)...);
      {%- endfor %}
      default : WALBERLA_ABORT("Shape type (" << ac.getShape(idx)->getShapeType() << ") could not be determined!");
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
