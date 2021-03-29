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

#include <mesa_pd/data/IContactAccessor.h>
#include <mesa_pd/data/ContactStorage.h>

#include <core/UniqueID.h>

#include <limits>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * @brief Basic ContactAccessor for the ContactStorage
 *
 * Provides get, set and getRef for all members of the ContactStorage.
 * Can be used as a basis class for a more advanced ContactAccessor.
 */
class ContactAccessor : public IContactAccessor
{
public:
   ContactAccessor(const std::shared_ptr<data::ContactStorage>& ps) : ps_(ps) {}
   ~ContactAccessor() override = default;

   {%- for prop in properties %}
   const {{prop.type}}& get{{prop.name | capFirst}}(const size_t p_idx) const {return ps_->get{{prop.name | capFirst}}(p_idx);}
   {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t p_idx) {return ps_->get{{prop.name | capFirst}}Ref(p_idx);}
   void set{{prop.name | capFirst}}(const size_t p_idx, const {{prop.type}}& v) { ps_->set{{prop.name | capFirst}}(p_idx, v);}
   {% endfor %}

   id_t getInvalidUid() const {return UniqueID<data::Contact>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of Contact specified by uid.
   * @param uid unique id of the Contact to be looked up
   * @return the index of the Contact or std::numeric_limits<size_t>::max() if the Contact is not found
   */
   size_t uidToIdx(const id_t& uid) const {auto it = ps_->find(uid); return it != ps_->end() ? it.getIdx() : std::numeric_limits<size_t>::max();}
   size_t size() const { return ps_->size(); }

   inline size_t create(const id_t& uid);
   inline size_t erase(const size_t& idx);
   inline size_t find(const id_t& uid);
protected:
   std::shared_ptr<data::ContactStorage> ps_;
};

inline size_t ContactAccessor::create(const id_t& uid)
{
   auto it = ps_->create(uid);
   return it.getIdx();
}
inline size_t ContactAccessor::erase(const size_t& idx)
{
   data::ContactStorage::iterator it(ps_.get(), idx);
   it = ps_->erase(it);
   return it.getIdx();
}
inline size_t ContactAccessor::find(const id_t& uid)
{
   auto it = ps_->find(uid);
   return it.getIdx();
}

/**
 * @brief Basic ContactAccessor which emulates a single Contact in a ContactStorage
 * without actually needing a ContactStorage. This class is used mainly for testing purposes.
 *
 * Provides get, set and getRef.
 */
class SingleContactAccessor : public IContactAccessor
{
public:
   ~SingleContactAccessor() override = default;

   {%- for prop in properties %}
   const {{prop.type}}& get{{prop.name | capFirst}}(const size_t /*p_idx*/) const {return {{prop.name}}_;}
   void set{{prop.name | capFirst}}(const size_t /*p_idx*/, const {{prop.type}}& v) { {{prop.name}}_ = v;}
   {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t /*p_idx*/) {return {{prop.name}}_;}
   {% endfor %}

   id_t getInvalidUid() const {return UniqueID<data::Contact>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of Contact specified by uid.
   * @param uid unique id of the Contact to be looked up
   * @return the index of the Contact or std::numeric_limits<size_t>::max() if the Contact is not found
   */
   size_t uidToIdx(const id_t& uid) const {return uid == uid_ ? 0 : std::numeric_limits<size_t>::max();}
   size_t size() const { return 1; }
private:
   {%- for prop in properties %}
   {{prop.type}} {{prop.name}}_;
   {%- endfor %}
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla
