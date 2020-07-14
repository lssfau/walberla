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
//! \file ParticleStorage.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <atomic>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <mesa_pd/data/ContactHistory.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>
{%- for include in includes %}
#include <{{include}}>
{%- endfor %}
#include <mesa_pd/data/STLOverloads.h>

#include <core/Abort.h>
#include <core/debug/Debug.h>
#include <core/math/AABB.h>
#include <core/mpi/MPIWrapper.h>
#include <core/OpenMP.h>
#include <core/STLIO.h>
#include <core/UniqueID.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class ParticleStorage;

class ParticleStorage
{
public:
   class Particle
   {
   public:
      constexpr Particle(ParticleStorage& storage, const size_t i) : storage_(storage), i_(i) {}
      constexpr Particle(const Particle&)  = default;
      constexpr Particle(Particle&&)  = default;

      Particle& operator=(const Particle& rhs);
      Particle& operator=(Particle&& rhs);

      Particle* operator->(){return this;}

      {% for prop in properties %}
      using {{prop.name}}_type = {{prop.type}};
      {%- endfor %}

      {% for prop in properties %}
      {{prop.name}}_type const & get{{prop.name | capFirst}}() const {return storage_.get{{prop.name | capFirst}}(i_);}
      {{prop.name}}_type& get{{prop.name | capFirst}}Ref() {return storage_.get{{prop.name | capFirst}}Ref(i_);}
      void set{{prop.name | capFirst}}({{prop.name}}_type const & v) { storage_.set{{prop.name | capFirst}}(i_, v);}
      {% endfor %}

      size_t getIdx() const {return i_;}
   public:
      ParticleStorage& storage_;
      const size_t i_;
   };

   class iterator
   {
   public:
      using iterator_category = std::random_access_iterator_tag;
      using value_type        = Particle;
      using pointer           = Particle*;
      using reference         = Particle&;
      using difference_type   = std::ptrdiff_t;

      explicit iterator(ParticleStorage* storage, const size_t i) : storage_(storage), i_(i) {}
      iterator(const iterator& it)         = default;
      iterator(iterator&& it)              = default;
      iterator& operator=(const iterator&) = default;
      iterator& operator=(iterator&&)      = default;


      Particle operator*() const {return Particle{*storage_, i_};}
      Particle operator->() const {return Particle{*storage_, i_};}
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
      ParticleStorage* storage_;
      size_t i_;
   };

   explicit ParticleStorage(const size_t size);

   iterator begin() { return iterator(this, 0); }
   iterator end()   { return iterator(this, size()); }
   Particle operator[](const size_t n) { return *iterator(this, n); }

   {% for prop in properties %}
   using {{prop.name}}_type = {{prop.type}};
   {%- endfor %}

   {% for prop in properties %}
   {{prop.name}}_type const & get{{prop.name | capFirst}}(const size_t idx) const {return {{prop.name}}_[idx];}
   {{prop.name}}_type& get{{prop.name | capFirst}}Ref(const size_t idx) {return {{prop.name}}_[idx];}
   void set{{prop.name | capFirst}}(const size_t idx, {{prop.name}}_type const & v) { {{prop.name}}_[idx] = v; }
   {% endfor %}

   /**
    * @brief creates a new particle and returns an iterator pointing to it
    *
    * \attention Use this function only if you know what you are doing!
    * Messing with the uid might break the simulation!
    * If you are unsure use create(bool) instead.
    * @param uid unique id of the particle to be created
    * @return iterator to the newly created particle
    */
   inline iterator create(const uid_type& uid);
   inline iterator create(const bool global = false);
   inline iterator erase(iterator& it);
   /// Finds the entry corresponding to \p uid.
   /// \return iterator to the object or end iterator
   inline iterator find(const uid_type& uid);
   inline void reserve(const size_t size);
   inline void clear();
   inline size_t size() const;
   template <class Compare>
   void sort(Compare comp);

   /**
    * Calls the provided functor \p func for all Particles.
    *
    * Additional arguments can be provided.
    * Call syntax for the provided functor
    * \code
    * func( *this, i, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticle(const bool openmp,
                               const Selector& selector,
                               Accessor& acForPS,
                               Func&& func,
                               Args&&... args);
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticle(const bool openmp,
                               const Selector& selector,
                               Accessor& acForPS,
                               Func&& func,
                               Args&&... args) const;
   /**
    * Calls the provided functor \p func for all Particle pairs.
    *
    * Additional arguments can be provided. No pairs with twice the same Particle.
    * Call syntax for the provided functor
    * \code
    * func( *this, i, j, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticlePair(const bool openmp,
                                   const Selector& selector,
                                   Accessor& acForPS,
                                   Func&& func,
                                   Args&&... args);
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticlePair(const bool openmp,
                                   const Selector& selector,
                                   Accessor& acForPS,
                                   Func&& func,
                                   Args&&... args) const;
   /**
    * Calls the provided functor \p func for all Particle pairs.
    *
    * Additional arguments can be provided. No pairs with twice the same Particle.
    * Index of the first particle i will be always smaller than j! No pair is called twice!
    * Call syntax for the provided functor
    * \code
    * func( *this, i, j, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticlePairHalf(const bool openmp,
                                       const Selector& selector,
                                       Accessor& acForPS,
                                       Func&& func,
                                       Args&&... args);
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   inline void forEachParticlePairHalf(const bool openmp,
                                       const Selector& selector,
                                       Accessor& acForPS,
                                       Func&& func,
                                       Args&&... args) const;

   private:
   {%- for prop in properties %}
   std::vector<{{prop.name}}_type> {{prop.name}}_ {};
   {%- endfor %}
   std::unordered_map<uid_type, size_t> uidToIdx_;
   static_assert(std::is_same<uid_type, id_t>::value,
                 "Property uid of type id_t is missing. This property is required!");
};
using Particle = ParticleStorage::Particle;

inline
ParticleStorage::Particle& ParticleStorage::Particle::operator=(const ParticleStorage::Particle& rhs)
{
   {%- for prop in properties %}
   get{{prop.name | capFirst}}Ref() = rhs.get{{prop.name | capFirst}}();
   {%- endfor %}
   return *this;
}

inline
ParticleStorage::Particle& ParticleStorage::Particle::operator=(ParticleStorage::Particle&& rhs)
{
   {%- for prop in properties %}
   get{{prop.name | capFirst}}Ref() = std::move(rhs.get{{prop.name | capFirst}}Ref());
   {%- endfor %}
   return *this;
}

inline
void swap(ParticleStorage::Particle lhs, ParticleStorage::Particle rhs)
{
   if (lhs.i_ == rhs.i_) return;
   {%- for prop in properties %}
   std::swap(lhs.get{{prop.name | capFirst}}Ref(), rhs.get{{prop.name | capFirst}}Ref());
   {%- endfor %}
}

inline
std::ostream& operator<<( std::ostream& os, const ParticleStorage::Particle& p )
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
ParticleStorage::iterator operator+(const ParticleStorage::iterator& it, const size_t n)
{
   return ParticleStorage::iterator(it.storage_, it.i_+n);
}

inline
ParticleStorage::iterator operator+(const size_t n, const ParticleStorage::iterator& it)
{
   return it + n;
}

inline
ParticleStorage::iterator operator-(const ParticleStorage::iterator& it, const size_t n)
{
   return ParticleStorage::iterator(it.storage_, it.i_-n);
}

inline
ParticleStorage::iterator::difference_type operator-(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   return int64_c(lhs.i_) - int64_c(rhs.i_);
}

inline bool operator==(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ == rhs.i_);
}
inline bool operator!=(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ != rhs.i_);
}
inline bool operator<(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ < rhs.i_);
}
inline bool operator>(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ > rhs.i_);
}
inline bool operator<=(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ <= rhs.i_);
}
inline bool operator>=(const ParticleStorage::iterator& lhs, const ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   return (lhs.i_ >= rhs.i_);
}

inline void swap(ParticleStorage::iterator& lhs, ParticleStorage::iterator& rhs)
{
   WALBERLA_ASSERT_EQUAL(lhs.storage_, rhs.storage_);
   std::swap(lhs.i_, rhs.i_);
}

inline
ParticleStorage::ParticleStorage(const size_t size)
{
   reserve(size);
}


inline ParticleStorage::iterator ParticleStorage::create(const id_t& uid)
{
   WALBERLA_ASSERT_EQUAL(uidToIdx_.find(uid),
                         uidToIdx_.end(),
                         "particle with the same uid(" << uid <<") already existing at index(" << uidToIdx_.find(uid)->second << ")");
   {%- for prop in properties %}
   {{prop.name}}_.emplace_back({{prop.defValue}});
   {%- endfor %}
   uid_.back() = uid;
   uidToIdx_[uid] = uid_.size() - 1;
   return iterator(this, size() - 1);
}

inline ParticleStorage::iterator ParticleStorage::create(const bool global)
{
   if (global)
   {
      auto it = create(UniqueID<Particle>::createGlobal());
      data::particle_flags::set(flags_.back(), data::particle_flags::GLOBAL);
      data::particle_flags::set(flags_.back(), data::particle_flags::NON_COMMUNICATING);
      return it;
   } else
   {
      return create(UniqueID<Particle>::create());
   }
}

inline ParticleStorage::iterator ParticleStorage::erase(iterator& it)
{
   //swap with last element and pop
   auto last = --end();
   auto numElementsRemoved = uidToIdx_.erase(it->getUid());
   WALBERLA_CHECK_EQUAL(numElementsRemoved,
                        1,
                        "Particle with uid " << it->getUid() << " cannot be removed (not existing).");
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

inline ParticleStorage::iterator ParticleStorage::find(const id_t& uid)
{
   //linear search through uid vector
   //auto it = std::find(uid_.begin(), uid_.end(), uid);
   //if (it == uid_.end()) return end();
   //return iterator(this, uint_c(std::distance(uid_.begin(), it)));

   //use unordered_map for faster lookup
   auto it = uidToIdx_.find(uid);
   if (it == uidToIdx_.end()) return end();
   WALBERLA_ASSERT_EQUAL(getUid(it->second), uid, "Lookup via uidToIdx map is not up to date!!!");
   return iterator(this, it->second);
}

inline void ParticleStorage::reserve(const size_t size)
{
   {%- for prop in properties %}
   {{prop.name}}_.reserve(size);
   {%- endfor %}
}

inline void ParticleStorage::clear()
{
   {%- for prop in properties %}
   {{prop.name}}_.clear();
   {%- endfor %}
   uidToIdx_.clear();
}

inline size_t ParticleStorage::size() const
{
   {%- for prop in properties %}
   //WALBERLA_ASSERT_EQUAL( {{properties[0].name}}_.size(), {{prop.name}}.size() );
   {%- endfor %}
   return {{properties[0].name}}_.size();
}

template <class Compare>
void ParticleStorage::sort(Compare comp)
{
   using WeightPair = std::pair<size_t, double>; //idx, weight

   const size_t length = size();
   std::vector<size_t>     newIdx(length); //where is old idx now?
   std::vector<size_t>     oldIdx(length); //what old idx is at that idx?
   std::vector<WeightPair> weight(length);

   for (size_t idx = 0; idx < length; ++idx)
   {
      newIdx[idx] = idx;
      oldIdx[idx] = idx;
      weight[idx] = std::make_pair(idx, comp.getWeight(operator[](idx)));
   }
   std::sort(weight.begin(), weight.end(), [](const WeightPair& lhs, const WeightPair& rhs){return lhs.second < rhs.second;});
   for (size_t idx = 0; idx < length; ++idx)
   {
      using std::swap;
      WALBERLA_ASSERT_IDENTICAL(weight[idx].second, comp.getWeight(operator[](newIdx[weight[idx].first])));

      WALBERLA_ASSERT_LESS_EQUAL(comp.getWeight(operator[](newIdx[weight[idx].first])), comp.getWeight(operator[](idx)));
      swap( operator[](idx), operator[](newIdx[weight[idx].first]) );
      const auto lhsIDX = idx;
      const auto rhsIDX = newIdx[weight[idx].first];

      newIdx[oldIdx[lhsIDX]] = rhsIDX;
      newIdx[oldIdx[rhsIDX]] = lhsIDX;

      swap(oldIdx[lhsIDX], oldIdx[rhsIDX]);
   }

   //rebuild lookup table
   uidToIdx_.clear();
   for (size_t idx = 0; idx < length; ++idx)
   {
      uidToIdx_[getUid(idx)] = idx;
   }
}

{%- for const in ["", "const"] %}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticle(const bool openmp,
                                             const Selector& selector,
                                             Accessor& acForPS,
                                             Func&& func, Args&&... args) {{const}}
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   {%- if module.enableOpenMP %}
   #pragma omp parallel for schedule(static) firstprivate(selector,func) if (openmp)
   {%- endif %}
   for (int64_t i = 0; i < int64_c(len); ++i)
   {
      if (selector(uint64_c(i), acForPS))
         func( uint64_c(i), std::forward<Args>(args)... );
   }
}
{%- endfor %}

{%- for const in ["", "const"] %}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePair(const bool openmp,
                                                 const Selector& selector,
                                                 Accessor& acForPS,
                                                 Func&& func, Args&&... args) {{const}}
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   {%- if module.enableOpenMP %}
   #pragma omp parallel for schedule(static) firstprivate(selector,func) if (openmp)
   {%- endif %}
   for (int64_t i = 0; i < int64_c(len); ++i)
   {
      for (int64_t j = 0; j < int64_c(len); ++j)
      {
         if (i!=j)
         {
            if (selector(uint64_c(i), uint64_c(j), acForPS))
            {
               func( uint64_c(i), uint64_c(j), std::forward<Args>(args)... );
            }
         }
      }
   }
}
{%- endfor %}

{%- for const in ["", "const"] %}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePairHalf(const bool openmp,
                                                     const Selector& selector,
                                                     Accessor& acForPS,
                                                     Func&& func, Args&&... args) {{const}}
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   {%- if module.enableOpenMP %}
   #pragma omp parallel for schedule(static) firstprivate(selector,func) if (openmp)
   {%- endif %}
   for (int64_t i = 0; i < int64_c(len); ++i)
   {
      for (int64_t j = i+1; j < int64_c(len); ++j)
      {
         if (selector(uint64_c(i), uint64_c(j), acForPS))
         {
            func( uint64_c(i), uint64_c(j), std::forward<Args>(args)... );
         }
      }
   }
}
{%- endfor %}


{%- for prop in properties %}
///Predicate that selects a certain property from a Particle
class SelectParticle{{prop.name | capFirst}}
{
public:
   using return_type = {{prop.type}};
   {{prop.type}}& operator()(data::Particle& p) const {return p.get{{prop.name | capFirst}}Ref();}
   {{prop.type}}& operator()(data::Particle&& p) const {return p.get{{prop.name | capFirst}}Ref();}
   {{prop.type}} const & operator()(const data::Particle& p) const {return p.get{{prop.name | capFirst}}();}
};
{%- endfor %}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
