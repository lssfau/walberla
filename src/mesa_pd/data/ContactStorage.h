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

      
      const walberla::id_t& getUid() const {return storage_.getUid(i_);}
      walberla::id_t& getUidRef() {return storage_.getUidRef(i_);}
      void setUid(const walberla::id_t& v) { storage_.setUid(i_, v);}
      
      const walberla::id_t& getId1() const {return storage_.getId1(i_);}
      walberla::id_t& getId1Ref() {return storage_.getId1Ref(i_);}
      void setId1(const walberla::id_t& v) { storage_.setId1(i_, v);}
      
      const walberla::id_t& getId2() const {return storage_.getId2(i_);}
      walberla::id_t& getId2Ref() {return storage_.getId2Ref(i_);}
      void setId2(const walberla::id_t& v) { storage_.setId2(i_, v);}
      
      const real_t& getDistance() const {return storage_.getDistance(i_);}
      real_t& getDistanceRef() {return storage_.getDistanceRef(i_);}
      void setDistance(const real_t& v) { storage_.setDistance(i_, v);}
      
      const walberla::mesa_pd::Vec3& getNormal() const {return storage_.getNormal(i_);}
      walberla::mesa_pd::Vec3& getNormalRef() {return storage_.getNormalRef(i_);}
      void setNormal(const walberla::mesa_pd::Vec3& v) { storage_.setNormal(i_, v);}
      
      const walberla::mesa_pd::Vec3& getPosition() const {return storage_.getPosition(i_);}
      walberla::mesa_pd::Vec3& getPositionRef() {return storage_.getPositionRef(i_);}
      void setPosition(const walberla::mesa_pd::Vec3& v) { storage_.setPosition(i_, v);}
      
      const walberla::mesa_pd::Vec3& getT() const {return storage_.getT(i_);}
      walberla::mesa_pd::Vec3& getTRef() {return storage_.getTRef(i_);}
      void setT(const walberla::mesa_pd::Vec3& v) { storage_.setT(i_, v);}
      
      const walberla::mesa_pd::Vec3& getO() const {return storage_.getO(i_);}
      walberla::mesa_pd::Vec3& getORef() {return storage_.getORef(i_);}
      void setO(const walberla::mesa_pd::Vec3& v) { storage_.setO(i_, v);}
      
      const walberla::mesa_pd::Vec3& getR1() const {return storage_.getR1(i_);}
      walberla::mesa_pd::Vec3& getR1Ref() {return storage_.getR1Ref(i_);}
      void setR1(const walberla::mesa_pd::Vec3& v) { storage_.setR1(i_, v);}
      
      const walberla::mesa_pd::Vec3& getR2() const {return storage_.getR2(i_);}
      walberla::mesa_pd::Vec3& getR2Ref() {return storage_.getR2Ref(i_);}
      void setR2(const walberla::mesa_pd::Vec3& v) { storage_.setR2(i_, v);}
      
      const real_t& getMu() const {return storage_.getMu(i_);}
      real_t& getMuRef() {return storage_.getMuRef(i_);}
      void setMu(const real_t& v) { storage_.setMu(i_, v);}
      
      const walberla::mesa_pd::Vec3& getP() const {return storage_.getP(i_);}
      walberla::mesa_pd::Vec3& getPRef() {return storage_.getPRef(i_);}
      void setP(const walberla::mesa_pd::Vec3& v) { storage_.setP(i_, v);}
      
      const walberla::mesa_pd::Mat3& getDiag_nto() const {return storage_.getDiag_nto(i_);}
      walberla::mesa_pd::Mat3& getDiag_ntoRef() {return storage_.getDiag_ntoRef(i_);}
      void setDiag_nto(const walberla::mesa_pd::Mat3& v) { storage_.setDiag_nto(i_, v);}
      
      const walberla::mesa_pd::Mat3& getDiag_nto_inv() const {return storage_.getDiag_nto_inv(i_);}
      walberla::mesa_pd::Mat3& getDiag_nto_invRef() {return storage_.getDiag_nto_invRef(i_);}
      void setDiag_nto_inv(const walberla::mesa_pd::Mat3& v) { storage_.setDiag_nto_inv(i_, v);}
      
      const walberla::mesa_pd::Mat2& getDiag_to_inv() const {return storage_.getDiag_to_inv(i_);}
      walberla::mesa_pd::Mat2& getDiag_to_invRef() {return storage_.getDiag_to_invRef(i_);}
      void setDiag_to_inv(const walberla::mesa_pd::Mat2& v) { storage_.setDiag_to_inv(i_, v);}
      
      const real_t& getDiag_n_inv() const {return storage_.getDiag_n_inv(i_);}
      real_t& getDiag_n_invRef() {return storage_.getDiag_n_invRef(i_);}
      void setDiag_n_inv(const real_t& v) { storage_.setDiag_n_inv(i_, v);}
      

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

   
   const walberla::id_t& getUid(const size_t idx) const {return uid_[idx];}
   walberla::id_t& getUidRef(const size_t idx) {return uid_[idx];}
   void setUid(const size_t idx, const walberla::id_t& v) { uid_[idx] = v; }
   
   const walberla::id_t& getId1(const size_t idx) const {return id1_[idx];}
   walberla::id_t& getId1Ref(const size_t idx) {return id1_[idx];}
   void setId1(const size_t idx, const walberla::id_t& v) { id1_[idx] = v; }
   
   const walberla::id_t& getId2(const size_t idx) const {return id2_[idx];}
   walberla::id_t& getId2Ref(const size_t idx) {return id2_[idx];}
   void setId2(const size_t idx, const walberla::id_t& v) { id2_[idx] = v; }
   
   const real_t& getDistance(const size_t idx) const {return distance_[idx];}
   real_t& getDistanceRef(const size_t idx) {return distance_[idx];}
   void setDistance(const size_t idx, const real_t& v) { distance_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getNormal(const size_t idx) const {return normal_[idx];}
   walberla::mesa_pd::Vec3& getNormalRef(const size_t idx) {return normal_[idx];}
   void setNormal(const size_t idx, const walberla::mesa_pd::Vec3& v) { normal_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t idx) const {return position_[idx];}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t idx) {return position_[idx];}
   void setPosition(const size_t idx, const walberla::mesa_pd::Vec3& v) { position_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getT(const size_t idx) const {return t_[idx];}
   walberla::mesa_pd::Vec3& getTRef(const size_t idx) {return t_[idx];}
   void setT(const size_t idx, const walberla::mesa_pd::Vec3& v) { t_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getO(const size_t idx) const {return o_[idx];}
   walberla::mesa_pd::Vec3& getORef(const size_t idx) {return o_[idx];}
   void setO(const size_t idx, const walberla::mesa_pd::Vec3& v) { o_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getR1(const size_t idx) const {return r1_[idx];}
   walberla::mesa_pd::Vec3& getR1Ref(const size_t idx) {return r1_[idx];}
   void setR1(const size_t idx, const walberla::mesa_pd::Vec3& v) { r1_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getR2(const size_t idx) const {return r2_[idx];}
   walberla::mesa_pd::Vec3& getR2Ref(const size_t idx) {return r2_[idx];}
   void setR2(const size_t idx, const walberla::mesa_pd::Vec3& v) { r2_[idx] = v; }
   
   const real_t& getMu(const size_t idx) const {return mu_[idx];}
   real_t& getMuRef(const size_t idx) {return mu_[idx];}
   void setMu(const size_t idx, const real_t& v) { mu_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getP(const size_t idx) const {return p_[idx];}
   walberla::mesa_pd::Vec3& getPRef(const size_t idx) {return p_[idx];}
   void setP(const size_t idx, const walberla::mesa_pd::Vec3& v) { p_[idx] = v; }
   
   const walberla::mesa_pd::Mat3& getDiag_nto(const size_t idx) const {return diag_nto_[idx];}
   walberla::mesa_pd::Mat3& getDiag_ntoRef(const size_t idx) {return diag_nto_[idx];}
   void setDiag_nto(const size_t idx, const walberla::mesa_pd::Mat3& v) { diag_nto_[idx] = v; }
   
   const walberla::mesa_pd::Mat3& getDiag_nto_inv(const size_t idx) const {return diag_nto_inv_[idx];}
   walberla::mesa_pd::Mat3& getDiag_nto_invRef(const size_t idx) {return diag_nto_inv_[idx];}
   void setDiag_nto_inv(const size_t idx, const walberla::mesa_pd::Mat3& v) { diag_nto_inv_[idx] = v; }
   
   const walberla::mesa_pd::Mat2& getDiag_to_inv(const size_t idx) const {return diag_to_inv_[idx];}
   walberla::mesa_pd::Mat2& getDiag_to_invRef(const size_t idx) {return diag_to_inv_[idx];}
   void setDiag_to_inv(const size_t idx, const walberla::mesa_pd::Mat2& v) { diag_to_inv_[idx] = v; }
   
   const real_t& getDiag_n_inv(const size_t idx) const {return diag_n_inv_[idx];}
   real_t& getDiag_n_invRef(const size_t idx) {return diag_n_inv_[idx];}
   void setDiag_n_inv(const size_t idx, const real_t& v) { diag_n_inv_[idx] = v; }
   

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
   std::vector<walberla::id_t> uid_ {};
   std::vector<walberla::id_t> id1_ {};
   std::vector<walberla::id_t> id2_ {};
   std::vector<real_t> distance_ {};
   std::vector<walberla::mesa_pd::Vec3> normal_ {};
   std::vector<walberla::mesa_pd::Vec3> position_ {};
   std::vector<walberla::mesa_pd::Vec3> t_ {};
   std::vector<walberla::mesa_pd::Vec3> o_ {};
   std::vector<walberla::mesa_pd::Vec3> r1_ {};
   std::vector<walberla::mesa_pd::Vec3> r2_ {};
   std::vector<real_t> mu_ {};
   std::vector<walberla::mesa_pd::Vec3> p_ {};
   std::vector<walberla::mesa_pd::Mat3> diag_nto_ {};
   std::vector<walberla::mesa_pd::Mat3> diag_nto_inv_ {};
   std::vector<walberla::mesa_pd::Mat2> diag_to_inv_ {};
   std::vector<real_t> diag_n_inv_ {};
   std::unordered_map<id_t, size_t> uidToIdx_;
   static_assert(std::is_same<decltype(uid_)::value_type, id_t>::value,
                 "Property uid of type id_t is missing. This property is required!");
};
using Contact = ContactStorage::Contact;

inline
ContactStorage::Contact& ContactStorage::Contact::operator=(const ContactStorage::Contact& rhs)
{
   getUidRef() = rhs.getUid();
   getId1Ref() = rhs.getId1();
   getId2Ref() = rhs.getId2();
   getDistanceRef() = rhs.getDistance();
   getNormalRef() = rhs.getNormal();
   getPositionRef() = rhs.getPosition();
   getTRef() = rhs.getT();
   getORef() = rhs.getO();
   getR1Ref() = rhs.getR1();
   getR2Ref() = rhs.getR2();
   getMuRef() = rhs.getMu();
   getPRef() = rhs.getP();
   getDiag_ntoRef() = rhs.getDiag_nto();
   getDiag_nto_invRef() = rhs.getDiag_nto_inv();
   getDiag_to_invRef() = rhs.getDiag_to_inv();
   getDiag_n_invRef() = rhs.getDiag_n_inv();
   return *this;
}

inline
ContactStorage::Contact& ContactStorage::Contact::operator=(ContactStorage::Contact&& rhs)
{
   getUidRef() = std::move(rhs.getUidRef());
   getId1Ref() = std::move(rhs.getId1Ref());
   getId2Ref() = std::move(rhs.getId2Ref());
   getDistanceRef() = std::move(rhs.getDistanceRef());
   getNormalRef() = std::move(rhs.getNormalRef());
   getPositionRef() = std::move(rhs.getPositionRef());
   getTRef() = std::move(rhs.getTRef());
   getORef() = std::move(rhs.getORef());
   getR1Ref() = std::move(rhs.getR1Ref());
   getR2Ref() = std::move(rhs.getR2Ref());
   getMuRef() = std::move(rhs.getMuRef());
   getPRef() = std::move(rhs.getPRef());
   getDiag_ntoRef() = std::move(rhs.getDiag_ntoRef());
   getDiag_nto_invRef() = std::move(rhs.getDiag_nto_invRef());
   getDiag_to_invRef() = std::move(rhs.getDiag_to_invRef());
   getDiag_n_invRef() = std::move(rhs.getDiag_n_invRef());
   return *this;
}

inline
std::ostream& operator<<( std::ostream& os, const ContactStorage::Contact& p )
{
   os << "==========    ==========" << "\n" <<
         "idx                 : " << p.getIdx() << "\n" <<
         "uid                 : " << p.getUid() << "\n" <<
         "id1                 : " << p.getId1() << "\n" <<
         "id2                 : " << p.getId2() << "\n" <<
         "distance            : " << p.getDistance() << "\n" <<
         "normal              : " << p.getNormal() << "\n" <<
         "position            : " << p.getPosition() << "\n" <<
         "t                   : " << p.getT() << "\n" <<
         "o                   : " << p.getO() << "\n" <<
         "r1                  : " << p.getR1() << "\n" <<
         "r2                  : " << p.getR2() << "\n" <<
         "mu                  : " << p.getMu() << "\n" <<
         "p                   : " << p.getP() << "\n" <<
         "diag_nto            : " << p.getDiag_nto() << "\n" <<
         "diag_nto_inv        : " << p.getDiag_nto_inv() << "\n" <<
         "diag_to_inv         : " << p.getDiag_to_inv() << "\n" <<
         "diag_n_inv          : " << p.getDiag_n_inv() << "\n" <<
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
   uid_.emplace_back(walberla::id_t(-1));
   id1_.emplace_back(walberla::id_t(-1));
   id2_.emplace_back(walberla::id_t(-1));
   distance_.emplace_back(real_t(1));
   normal_.emplace_back(real_t(0));
   position_.emplace_back(real_t(0));
   t_.emplace_back(real_t(0));
   o_.emplace_back(real_t(0));
   r1_.emplace_back(real_t(0));
   r2_.emplace_back(real_t(0));
   mu_.emplace_back(real_t(0));
   p_.emplace_back(real_t(0));
   diag_nto_.emplace_back(real_t(0));
   diag_nto_inv_.emplace_back(real_t(0));
   diag_to_inv_.emplace_back(real_t(0));
   diag_n_inv_.emplace_back(real_t(0));
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
   uid_.pop_back();
   id1_.pop_back();
   id2_.pop_back();
   distance_.pop_back();
   normal_.pop_back();
   position_.pop_back();
   t_.pop_back();
   o_.pop_back();
   r1_.pop_back();
   r2_.pop_back();
   mu_.pop_back();
   p_.pop_back();
   diag_nto_.pop_back();
   diag_nto_inv_.pop_back();
   diag_to_inv_.pop_back();
   diag_n_inv_.pop_back();
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
   uid_.reserve(size);
   id1_.reserve(size);
   id2_.reserve(size);
   distance_.reserve(size);
   normal_.reserve(size);
   position_.reserve(size);
   t_.reserve(size);
   o_.reserve(size);
   r1_.reserve(size);
   r2_.reserve(size);
   mu_.reserve(size);
   p_.reserve(size);
   diag_nto_.reserve(size);
   diag_nto_inv_.reserve(size);
   diag_to_inv_.reserve(size);
   diag_n_inv_.reserve(size);
}

inline void ContactStorage::clear()
{
   uid_.clear();
   id1_.clear();
   id2_.clear();
   distance_.clear();
   normal_.clear();
   position_.clear();
   t_.clear();
   o_.clear();
   r1_.clear();
   r2_.clear();
   mu_.clear();
   p_.clear();
   diag_nto_.clear();
   diag_nto_inv_.clear();
   diag_to_inv_.clear();
   diag_n_inv_.clear();
   uidToIdx_.clear();
}

inline size_t ContactStorage::size() const
{
   //WALBERLA_ASSERT_EQUAL( uid_.size(), uid.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), id1.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), id2.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), distance.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), normal.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), position.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), t.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), o.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), r1.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), r2.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), mu.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), p.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), diag_nto.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), diag_nto_inv.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), diag_to_inv.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), diag_n_inv.size() );
   return uid_.size();
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ContactStorage::forEachContact(const bool openmp, const Selector& selector,
                                           Accessor& acForPS, Func&& func, Args&&... args) 
{
   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   for (int64_t i = 0; i < int64_c(len); ++i)
      if (selector(uint64_c(i), acForPS)){
         func( uint64_c(i), std::forward<Args>(args)... );
      }
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ContactStorage::forEachContact(const bool openmp, const Selector& selector,
                                           Accessor& acForPS, Func&& func, Args&&... args) const
{
   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   for (int64_t i = 0; i < int64_c(len); ++i)
      if (selector(uint64_c(i), acForPS)){
         func( uint64_c(i), std::forward<Args>(args)... );
      }
}
///Predicate that selects a certain property from a Contact
class SelectContactUid
{
public:
   using return_type = walberla::id_t;
   walberla::id_t& operator()(data::Contact& p) const {return p.getUidRef();}
   walberla::id_t& operator()(data::Contact&& p) const {return p.getUidRef();}
   const walberla::id_t& operator()(const data::Contact& p) const {return p.getUid();}
};
///Predicate that selects a certain property from a Contact
class SelectContactId1
{
public:
   using return_type = walberla::id_t;
   walberla::id_t& operator()(data::Contact& p) const {return p.getId1Ref();}
   walberla::id_t& operator()(data::Contact&& p) const {return p.getId1Ref();}
   const walberla::id_t& operator()(const data::Contact& p) const {return p.getId1();}
};
///Predicate that selects a certain property from a Contact
class SelectContactId2
{
public:
   using return_type = walberla::id_t;
   walberla::id_t& operator()(data::Contact& p) const {return p.getId2Ref();}
   walberla::id_t& operator()(data::Contact&& p) const {return p.getId2Ref();}
   const walberla::id_t& operator()(const data::Contact& p) const {return p.getId2();}
};
///Predicate that selects a certain property from a Contact
class SelectContactDistance
{
public:
   using return_type = real_t;
   real_t& operator()(data::Contact& p) const {return p.getDistanceRef();}
   real_t& operator()(data::Contact&& p) const {return p.getDistanceRef();}
   const real_t& operator()(const data::Contact& p) const {return p.getDistance();}
};
///Predicate that selects a certain property from a Contact
class SelectContactNormal
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getNormalRef();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getNormalRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getNormal();}
};
///Predicate that selects a certain property from a Contact
class SelectContactPosition
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getPositionRef();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getPositionRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getPosition();}
};
///Predicate that selects a certain property from a Contact
class SelectContactT
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getTRef();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getTRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getT();}
};
///Predicate that selects a certain property from a Contact
class SelectContactO
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getORef();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getORef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getO();}
};
///Predicate that selects a certain property from a Contact
class SelectContactR1
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getR1Ref();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getR1Ref();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getR1();}
};
///Predicate that selects a certain property from a Contact
class SelectContactR2
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getR2Ref();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getR2Ref();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getR2();}
};
///Predicate that selects a certain property from a Contact
class SelectContactMu
{
public:
   using return_type = real_t;
   real_t& operator()(data::Contact& p) const {return p.getMuRef();}
   real_t& operator()(data::Contact&& p) const {return p.getMuRef();}
   const real_t& operator()(const data::Contact& p) const {return p.getMu();}
};
///Predicate that selects a certain property from a Contact
class SelectContactP
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Contact& p) const {return p.getPRef();}
   walberla::mesa_pd::Vec3& operator()(data::Contact&& p) const {return p.getPRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Contact& p) const {return p.getP();}
};
///Predicate that selects a certain property from a Contact
class SelectContactDiag_nto
{
public:
   using return_type = walberla::mesa_pd::Mat3;
   walberla::mesa_pd::Mat3& operator()(data::Contact& p) const {return p.getDiag_ntoRef();}
   walberla::mesa_pd::Mat3& operator()(data::Contact&& p) const {return p.getDiag_ntoRef();}
   const walberla::mesa_pd::Mat3& operator()(const data::Contact& p) const {return p.getDiag_nto();}
};
///Predicate that selects a certain property from a Contact
class SelectContactDiag_nto_inv
{
public:
   using return_type = walberla::mesa_pd::Mat3;
   walberla::mesa_pd::Mat3& operator()(data::Contact& p) const {return p.getDiag_nto_invRef();}
   walberla::mesa_pd::Mat3& operator()(data::Contact&& p) const {return p.getDiag_nto_invRef();}
   const walberla::mesa_pd::Mat3& operator()(const data::Contact& p) const {return p.getDiag_nto_inv();}
};
///Predicate that selects a certain property from a Contact
class SelectContactDiag_to_inv
{
public:
   using return_type = walberla::mesa_pd::Mat2;
   walberla::mesa_pd::Mat2& operator()(data::Contact& p) const {return p.getDiag_to_invRef();}
   walberla::mesa_pd::Mat2& operator()(data::Contact&& p) const {return p.getDiag_to_invRef();}
   const walberla::mesa_pd::Mat2& operator()(const data::Contact& p) const {return p.getDiag_to_inv();}
};
///Predicate that selects a certain property from a Contact
class SelectContactDiag_n_inv
{
public:
   using return_type = real_t;
   real_t& operator()(data::Contact& p) const {return p.getDiag_n_invRef();}
   real_t& operator()(data::Contact&& p) const {return p.getDiag_n_invRef();}
   const real_t& operator()(const data::Contact& p) const {return p.getDiag_n_inv();}
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla