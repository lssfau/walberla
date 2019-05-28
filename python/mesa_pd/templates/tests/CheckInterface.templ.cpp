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
//! \file {{InterfaceTestName}}.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/{{KernelInclude}}>

#include <core/UniqueID.h>

#include <map>

namespace walberla {
namespace mesa_pd {

class Accessor : public data::IAccessor
{
public:
   virtual ~Accessor() = default;

   {%- for prop in interface %}
   {%- if 'g' in prop.access %}
   const {{prop.type}}& get{{prop.name | capFirst}}(const size_t /*p_idx*/) const {return {{prop.name}}_;}
   {%- endif %}
   {%- if 's' in prop.access %}
   void set{{prop.name | capFirst}}(const size_t /*p_idx*/, const {{prop.type}}& v) { {{prop.name}}_ = v;}
   {%- endif %}
   {%- if 'r' in prop.access %}
   {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t /*p_idx*/) {return {{prop.name}}_;}
   {%- endif %}
   {% endfor %}

   id_t getInvalidUid() const {return UniqueID<int>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of particle specified by uid.
   * @param uid unique id of the particle to be looked up
   * @return the index of the particle or std::numeric_limits<size_t>::max() if the particle is not found
   */
   size_t uidToIdx(const id_t& /*uid*/) const {return 0;}
   size_t size() const { return 1; }
private:
   {%- for prop in interface %}
   {{prop.type}} {{prop.name}}_;
   {%- endfor %}
};

{{ExplicitInstantiation}}

} //namespace mesa_pd
} //namespace walberla
