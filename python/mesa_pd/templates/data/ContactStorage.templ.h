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
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================

/**
 * Storage for detected contacts which can be used to perform actions
 * for all contacts, e.g. resolving, relaxation schemes...
 * The InsertIntoContactStorage-Kernel can be used to insert a contact.
 */

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
{%- for include in includes %}
#include <{{include}}>
{%- endfor %}
#include <mesa_pd/data/STLOverloads.h>

#include <core/Abort.h>
#include <core/debug/Debug.h>
#include <core/math/AABB.h>
#include <core/OpenMP.h>
#include <core/STLIO.h>
#include <core/UniqueID.h>

#include <atomic>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace data {

class ContactStorage;

class ContactStorage
{
public:
   class Contact
   {
   public:
      constexpr Contact(ContactStorage& storage, const size_t i) : storage_(storage), i_(i) {}
      constexpr Contact(const Contact&)  = default;
      constexpr Contact(Contact&&)  = default;

      Contact& operator=(const Contact& rhs);
      Contact& operator=(Contact&& rhs);

      Contact* operator->(){return this;}

      {% for prop in properties %}
      const {{prop.type}}& get{{prop.name | capFirst}}() const {return storage_.get{{prop.name | capFirst}}(i_);}
      {{prop.type}}& get{{prop.name | capFirst}}Ref() {return storage_.get{{prop.name | capFirst}}Ref(i_);}
      void set{{prop.name | capFirst}}(const {{prop.type}}& v) { storage_.set{{prop.name | capFirst}}(i_, v);}
      {% endfor %}

      size_t getIdx() const {return i_;}
   public:
      ContactStorage& storage_;
      const size_t i_;
   };

   class iterator
   {
   public:
      using iterator_category = std::random_access_iterator_tag;
      using value_type        = Contact;
      using pointer           = Contact*;
      using reference         = Contact&;
      using difference_type   = std::ptrdiff_t;

      explicit iterator(ContactStorage* storage, const size_t i) : storage_(storage), i_(i) {}
      iterator(const iterator& it)         = default;
      iterator(iterator&& it)              = default;
      iterator& operator=(const iterator&) = default;
      iterator& operator=(iterator&&)      = default;


      Contact operator*(){return Contact{*storage_, i_};}
      Contact operator->(){return Contact{*storage_, i_};}
      iterator& operator++(){ ++i_; return *this; }
      iterator operator++(int){ iterator tmp(*this); ++(*this); return tmp; }
      iterator& operator--(){ --i_; return *this; }
      iterator operator--(int){ iterator tmp(*this); --(*this); return tmp; }

      iterator& operator+=(const size_t n){ i_+=n; return *this; }
      iterator& operator-=(const size_t n){ i_-=n; return *this; }

      friend iterator operator+(const iterator& it, const size_t n);
      friend iterator operator+(const size_t n, const iterator& it);
      friend iterator operator-(const iterator& it, const size_t n);
      friend difference_type operator-(const iterator& lhs, const iterator& rhs);

      friend bool operator==(const iterator& lhs, const iterator& rhs);
      friend bool operator!=(const iterator& lhs, const iterator& rhs);
      friend bool operator<(const iterator& lhs, const iterator& rhs);
      friend bool operator>(const iterator& lhs, const iterator& rhs);
      friend bool operator<=(const iterator& lhs, const iterator& rhs);
      friend bool operator>=(const iterator& lhs, const iterator& rhs);

      friend void swap(iterator& lhs, iterator& rhs);

      size_t getIdx() const {return i_;}
   private:
      ContactStorage* storage_;
      size_t i_;
   };

   explicit ContactStorage(const size_t size);

   iterator begin() { return iterator(this, 0); }
   iterator end()   { return iterator(this, size()); }
   iterator operator[](const size_t n) { return iterator(this, n); }

   {% for prop in properties %}
   const {{prop.type}}& get{{prop.name | capFirst}}(const size_t idx) const {return {{prop.name}}_[idx];}
   {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t idx) {return {{prop.name}}_[idx];}
   void set{{prop.name | capFirst}}(const size_t idx, const {{prop.type}}& v) { {{prop.name}}_[idx] = v; }
   {% endfor %}

   /**
    * @brief creates a new Contact and returns an iterator pointing to it
    *
    * \attention Use this function only if you know what you are doing!
    * Messing with the uid might break the simulation!
    * If you are unsure use create(bool) instead.
    * @param uid unique id of the Contact to be created
    * @return iterator to the newly created Contact
    */
   inline iterator create(const id_t& uid);
   inline iterator create();
   inline iterator erase(iterator& it);
   /// Finds the entry corresponding to \p uid.
   /// \return iterator to the object or end iterator
   inline iterator find(const id_t& uid);
   inline void reserve(const size_t size);
   inline void clear();
   inline size_t size() const;

   /**
    * Calls the provided functor \p func for all Contacts selected by the selector.
    *
    * Additional arguments can be provided.
    * Call syntax for the provided functor
    * \code
    * func( *this, i, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachContact(const bool openmp, const Selector& selector,
                              Accessor& acForPS, Func&& func, Args&&... args);
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachContact(const bool openmp, const Selector& selector,
                              Accessor& acForPS, Func&& func, Args&&... args) const;

   private:
   {%- for prop in properties %}
   std::vector<{{prop.type}}> {{prop.name}}_ {};
   {%- endfor %}
   std::unordered_map<id_t, size_t> uidToIdx_;
   static_assert(std::is_same<decltype(uid_)::value_type, id_t>::value,
                 "Property uid of type id_t is missing. This property is required!");
};
using Contact = ContactStorage::Contact;

inline
ContactStorage::Contact& ContactStorage::Contact::operator=(const ContactStorage::Contact& rhs)
{
   {%- for prop in properties %}
   get{{prop.name | capFirst}}Ref() = rhs.get{{prop.name | capFirst}}();
   {%- endfor %}
   return *this;
}

inline
ContactStorage::Contact& ContactStorage::Contact::operator=(ContactStorage::Contact&& rhs)
{
   {%- for prop in properties %}
   get{{prop.name | capFirst}}Ref() = std::move(rhs.get{{prop.name | capFirst}}Ref());
   {%- endfor %}
   return *this;
}

inline
std::ostream& operator<<( std::ostream& os, const ContactStorage::Contact& p )
{
   os << "==========  {{StorageType | upper}}  ==========" << "\n" <<
         "idx                 : " << p.getIdx() << "\n" <<
   {%- for prop in properties %}
         "{{'%-20s'|format(prop.name)}}: " << p.get{{prop.name | capFirst}}() << "\n" <<
   {%- endfor %}
         "================================" << std::endl;
   return os;
}

inline
ContactStorage::iterator operator+(const ContactStorage::iterator& it, const size_t n)
{
   return ContactStorage::iterator(it.storage_, it.i_+n);
}

inline
ContactStorage::iterator operator+(const size_t n, const ContactStorage::iterator& it)
{
   return it + n;
}

inline
ContactStorage::iterator operator-(const ContactStorage::iterator& it, const size_t n)
{
   return ContactStorage::iterator(it.storage_, it.i_-n);
}

inline
ContactStorage::iterator::difference_type operator-(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   return int64_c(lhs.i_) - int64_c(rhs.i_);
}

inline bool operator==(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ == rhs.i_);
}
inline bool operator!=(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ != rhs.i_);
}
inline bool operator<(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ < rhs.i_);
}
inline bool operator>(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ > rhs.i_);
}
inline bool operator<=(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ <= rhs.i_);
}
inline bool operator>=(const ContactStorage::iterator& lhs, const ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ >= rhs.i_);
}

inline void swap(ContactStorage::iterator& lhs, ContactStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   std::swap(lhs.i_, rhs.i_);
}

inline
ContactStorage::ContactStorage(const size_t size)
{
   reserve(size);
}


inline ContactStorage::iterator ContactStorage::create(const id_t& uid)
{
   WALBERLA_ASSERT_EQUAL(uidToIdx_.find(uid),
                         uidToIdx_.end(),
                         "Contact with the same uid(" << uid <<") already existing at index(" << uidToIdx_.find(uid)->second << ")");
   {%- for prop in properties %}
   {{prop.name}}_.emplace_back({{prop.defValue}});
   {%- endfor %}
   uid_.back() = uid;
   uidToIdx_[uid] = uid_.size() - 1;
   return iterator(this, size() - 1);
}

inline ContactStorage::iterator ContactStorage::create()
{
   return create(UniqueID<Contact>::create());
}

inline ContactStorage::iterator ContactStorage::erase(iterator& it)
{
   //swap with last element and pop
   auto last = --end();
   auto numElementsRemoved = uidToIdx_.erase(it->getUid());
   WALBERLA_CHECK_EQUAL(numElementsRemoved,
                        1,
                        "Contact with uid " << it->getUid() << " cannot be removed (not existing).");
   if (it != last) //skip swap if last element is removed
   {
      *it = *last;
      uidToIdx_[it->getUid()] = it.getIdx();
   }
   {%- for prop in properties %}
   {{prop.name}}_.pop_back();
   {%- endfor %}
   return it;
}

inline ContactStorage::iterator ContactStorage::find(const id_t& uid)
{
   //linear search through uid vector
   //auto it = std::find(uid_.begin(), uid_.end(), uid);
   //if (it == uid_.end()) return end();
   //return iterator(this, uint_c(std::distance(uid_.begin(), it)));

   //use unordered_map for faster lookup
   auto it = uidToIdx_.find(uid);
   if (it == uidToIdx_.end()) return end();
   WALBERLA_ASSERT_EQUAL(it->first, uid, "Lookup via uidToIdx map is not up to date!!!");
   return iterator(this, it->second);
}

inline void ContactStorage::reserve(const size_t size)
{
   {%- for prop in properties %}
   {{prop.name}}_.reserve(size);
   {%- endfor %}
}

inline void ContactStorage::clear()
{
   {%- for prop in properties %}
   {{prop.name}}_.clear();
   {%- endfor %}
   uidToIdx_.clear();
}

inline size_t ContactStorage::size() const
{
   {%- for prop in properties %}
   //WALBERLA_ASSERT_EQUAL( {{properties[0].name}}_.size(), {{prop.name}}.size() );
   {%- endfor %}
   return {{properties[0].name}}_.size();
}

{%- for const in ["", "const"] %}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ContactStorage::forEachContact(const bool openmp, const Selector& selector,
                                           Accessor& acForPS, Func&& func, Args&&... args) {{const}}
{
   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   {%- if module.enableOpenMP %}
   #pragma omp parallel for schedule(static) firstprivate(selector,func) if (openmp)
   {%- endif %}
   for (int64_t i = 0; i < int64_c(len); ++i)
      if (selector(uint64_c(i), acForPS)){
         func( uint64_c(i), std::forward<Args>(args)... );
      }
}
{%- endfor %}


{%- for prop in properties %}
///Predicate that selects a certain property from a Contact
class SelectContact{{prop.name | capFirst}}
{
public:
   using return_type = {{prop.type}};
   {{prop.type}}& operator()(data::Contact& p) const {return p.get{{prop.name | capFirst}}Ref();}
   {{prop.type}}& operator()(data::Contact&& p) const {return p.get{{prop.name | capFirst}}Ref();}
   const {{prop.type}}& operator()(const data::Contact& p) const {return p.get{{prop.name | capFirst}}();}
};
{%- endfor %}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
