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
#include <vector>

#include <mesa_pd/data/ContactHistory.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/STLOverloads.h>

#include <core/Abort.h>
#include <core/debug/Debug.h>
#include <core/math/AABB.h>
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

      
      const walberla::id_t& getUid() const {return storage_.getUid(i_);}
      walberla::id_t& getUidRef() {return storage_.getUidRef(i_);}
      void setUid(const walberla::id_t& v) { storage_.setUid(i_, v);}
      
      const walberla::mesa_pd::Vec3& getPosition() const {return storage_.getPosition(i_);}
      walberla::mesa_pd::Vec3& getPositionRef() {return storage_.getPositionRef(i_);}
      void setPosition(const walberla::mesa_pd::Vec3& v) { storage_.setPosition(i_, v);}
      
      const walberla::real_t& getInteractionRadius() const {return storage_.getInteractionRadius(i_);}
      walberla::real_t& getInteractionRadiusRef() {return storage_.getInteractionRadiusRef(i_);}
      void setInteractionRadius(const walberla::real_t& v) { storage_.setInteractionRadius(i_, v);}
      
      const walberla::mesa_pd::data::particle_flags::FlagT& getFlags() const {return storage_.getFlags(i_);}
      walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef() {return storage_.getFlagsRef(i_);}
      void setFlags(const walberla::mesa_pd::data::particle_flags::FlagT& v) { storage_.setFlags(i_, v);}
      
      const int& getOwner() const {return storage_.getOwner(i_);}
      int& getOwnerRef() {return storage_.getOwnerRef(i_);}
      void setOwner(const int& v) { storage_.setOwner(i_, v);}
      
      const std::vector<int>& getGhostOwners() const {return storage_.getGhostOwners(i_);}
      std::vector<int>& getGhostOwnersRef() {return storage_.getGhostOwnersRef(i_);}
      void setGhostOwners(const std::vector<int>& v) { storage_.setGhostOwners(i_, v);}
      
      const size_t& getShapeID() const {return storage_.getShapeID(i_);}
      size_t& getShapeIDRef() {return storage_.getShapeIDRef(i_);}
      void setShapeID(const size_t& v) { storage_.setShapeID(i_, v);}
      
      const walberla::mesa_pd::Rot3& getRotation() const {return storage_.getRotation(i_);}
      walberla::mesa_pd::Rot3& getRotationRef() {return storage_.getRotationRef(i_);}
      void setRotation(const walberla::mesa_pd::Rot3& v) { storage_.setRotation(i_, v);}
      
      const walberla::mesa_pd::Vec3& getAngularVelocity() const {return storage_.getAngularVelocity(i_);}
      walberla::mesa_pd::Vec3& getAngularVelocityRef() {return storage_.getAngularVelocityRef(i_);}
      void setAngularVelocity(const walberla::mesa_pd::Vec3& v) { storage_.setAngularVelocity(i_, v);}
      
      const walberla::mesa_pd::Vec3& getTorque() const {return storage_.getTorque(i_);}
      walberla::mesa_pd::Vec3& getTorqueRef() {return storage_.getTorqueRef(i_);}
      void setTorque(const walberla::mesa_pd::Vec3& v) { storage_.setTorque(i_, v);}
      
      const walberla::mesa_pd::Vec3& getLinearVelocity() const {return storage_.getLinearVelocity(i_);}
      walberla::mesa_pd::Vec3& getLinearVelocityRef() {return storage_.getLinearVelocityRef(i_);}
      void setLinearVelocity(const walberla::mesa_pd::Vec3& v) { storage_.setLinearVelocity(i_, v);}
      
      const walberla::real_t& getInvMass() const {return storage_.getInvMass(i_);}
      walberla::real_t& getInvMassRef() {return storage_.getInvMassRef(i_);}
      void setInvMass(const walberla::real_t& v) { storage_.setInvMass(i_, v);}
      
      const walberla::mesa_pd::Vec3& getForce() const {return storage_.getForce(i_);}
      walberla::mesa_pd::Vec3& getForceRef() {return storage_.getForceRef(i_);}
      void setForce(const walberla::mesa_pd::Vec3& v) { storage_.setForce(i_, v);}
      
      const walberla::mesa_pd::Vec3& getOldForce() const {return storage_.getOldForce(i_);}
      walberla::mesa_pd::Vec3& getOldForceRef() {return storage_.getOldForceRef(i_);}
      void setOldForce(const walberla::mesa_pd::Vec3& v) { storage_.setOldForce(i_, v);}
      
      const walberla::mesa_pd::Vec3& getOldTorque() const {return storage_.getOldTorque(i_);}
      walberla::mesa_pd::Vec3& getOldTorqueRef() {return storage_.getOldTorqueRef(i_);}
      void setOldTorque(const walberla::mesa_pd::Vec3& v) { storage_.setOldTorque(i_, v);}
      
      const uint_t& getType() const {return storage_.getType(i_);}
      uint_t& getTypeRef() {return storage_.getTypeRef(i_);}
      void setType(const uint_t& v) { storage_.setType(i_, v);}
      
      const int& getNextParticle() const {return storage_.getNextParticle(i_);}
      int& getNextParticleRef() {return storage_.getNextParticleRef(i_);}
      void setNextParticle(const int& v) { storage_.setNextParticle(i_, v);}
      
      const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistory() const {return storage_.getOldContactHistory(i_);}
      std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef() {return storage_.getOldContactHistoryRef(i_);}
      void setOldContactHistory(const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { storage_.setOldContactHistory(i_, v);}
      
      const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistory() const {return storage_.getNewContactHistory(i_);}
      std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef() {return storage_.getNewContactHistoryRef(i_);}
      void setNewContactHistory(const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { storage_.setNewContactHistory(i_, v);}
      
      const walberla::real_t& getTemperature() const {return storage_.getTemperature(i_);}
      walberla::real_t& getTemperatureRef() {return storage_.getTemperatureRef(i_);}
      void setTemperature(const walberla::real_t& v) { storage_.setTemperature(i_, v);}
      
      const walberla::real_t& getHeatFlux() const {return storage_.getHeatFlux(i_);}
      walberla::real_t& getHeatFluxRef() {return storage_.getHeatFluxRef(i_);}
      void setHeatFlux(const walberla::real_t& v) { storage_.setHeatFlux(i_, v);}
      

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
   iterator operator[](const size_t n) { return iterator(this, n); }

   
   const walberla::id_t& getUid(const size_t idx) const {return uid_[idx];}
   walberla::id_t& getUidRef(const size_t idx) {return uid_[idx];}
   void setUid(const size_t idx, const walberla::id_t& v) { uid_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t idx) const {return position_[idx];}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t idx) {return position_[idx];}
   void setPosition(const size_t idx, const walberla::mesa_pd::Vec3& v) { position_[idx] = v; }
   
   const walberla::real_t& getInteractionRadius(const size_t idx) const {return interactionRadius_[idx];}
   walberla::real_t& getInteractionRadiusRef(const size_t idx) {return interactionRadius_[idx];}
   void setInteractionRadius(const size_t idx, const walberla::real_t& v) { interactionRadius_[idx] = v; }
   
   const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t idx) const {return flags_[idx];}
   walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef(const size_t idx) {return flags_[idx];}
   void setFlags(const size_t idx, const walberla::mesa_pd::data::particle_flags::FlagT& v) { flags_[idx] = v; }
   
   const int& getOwner(const size_t idx) const {return owner_[idx];}
   int& getOwnerRef(const size_t idx) {return owner_[idx];}
   void setOwner(const size_t idx, const int& v) { owner_[idx] = v; }
   
   const std::vector<int>& getGhostOwners(const size_t idx) const {return ghostOwners_[idx];}
   std::vector<int>& getGhostOwnersRef(const size_t idx) {return ghostOwners_[idx];}
   void setGhostOwners(const size_t idx, const std::vector<int>& v) { ghostOwners_[idx] = v; }
   
   const size_t& getShapeID(const size_t idx) const {return shapeID_[idx];}
   size_t& getShapeIDRef(const size_t idx) {return shapeID_[idx];}
   void setShapeID(const size_t idx, const size_t& v) { shapeID_[idx] = v; }
   
   const walberla::mesa_pd::Rot3& getRotation(const size_t idx) const {return rotation_[idx];}
   walberla::mesa_pd::Rot3& getRotationRef(const size_t idx) {return rotation_[idx];}
   void setRotation(const size_t idx, const walberla::mesa_pd::Rot3& v) { rotation_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t idx) const {return angularVelocity_[idx];}
   walberla::mesa_pd::Vec3& getAngularVelocityRef(const size_t idx) {return angularVelocity_[idx];}
   void setAngularVelocity(const size_t idx, const walberla::mesa_pd::Vec3& v) { angularVelocity_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getTorque(const size_t idx) const {return torque_[idx];}
   walberla::mesa_pd::Vec3& getTorqueRef(const size_t idx) {return torque_[idx];}
   void setTorque(const size_t idx, const walberla::mesa_pd::Vec3& v) { torque_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t idx) const {return linearVelocity_[idx];}
   walberla::mesa_pd::Vec3& getLinearVelocityRef(const size_t idx) {return linearVelocity_[idx];}
   void setLinearVelocity(const size_t idx, const walberla::mesa_pd::Vec3& v) { linearVelocity_[idx] = v; }
   
   const walberla::real_t& getInvMass(const size_t idx) const {return invMass_[idx];}
   walberla::real_t& getInvMassRef(const size_t idx) {return invMass_[idx];}
   void setInvMass(const size_t idx, const walberla::real_t& v) { invMass_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getForce(const size_t idx) const {return force_[idx];}
   walberla::mesa_pd::Vec3& getForceRef(const size_t idx) {return force_[idx];}
   void setForce(const size_t idx, const walberla::mesa_pd::Vec3& v) { force_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getOldForce(const size_t idx) const {return oldForce_[idx];}
   walberla::mesa_pd::Vec3& getOldForceRef(const size_t idx) {return oldForce_[idx];}
   void setOldForce(const size_t idx, const walberla::mesa_pd::Vec3& v) { oldForce_[idx] = v; }
   
   const walberla::mesa_pd::Vec3& getOldTorque(const size_t idx) const {return oldTorque_[idx];}
   walberla::mesa_pd::Vec3& getOldTorqueRef(const size_t idx) {return oldTorque_[idx];}
   void setOldTorque(const size_t idx, const walberla::mesa_pd::Vec3& v) { oldTorque_[idx] = v; }
   
   const uint_t& getType(const size_t idx) const {return type_[idx];}
   uint_t& getTypeRef(const size_t idx) {return type_[idx];}
   void setType(const size_t idx, const uint_t& v) { type_[idx] = v; }
   
   const int& getNextParticle(const size_t idx) const {return nextParticle_[idx];}
   int& getNextParticleRef(const size_t idx) {return nextParticle_[idx];}
   void setNextParticle(const size_t idx, const int& v) { nextParticle_[idx] = v; }
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistory(const size_t idx) const {return oldContactHistory_[idx];}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef(const size_t idx) {return oldContactHistory_[idx];}
   void setOldContactHistory(const size_t idx, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { oldContactHistory_[idx] = v; }
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistory(const size_t idx) const {return newContactHistory_[idx];}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef(const size_t idx) {return newContactHistory_[idx];}
   void setNewContactHistory(const size_t idx, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { newContactHistory_[idx] = v; }
   
   const walberla::real_t& getTemperature(const size_t idx) const {return temperature_[idx];}
   walberla::real_t& getTemperatureRef(const size_t idx) {return temperature_[idx];}
   void setTemperature(const size_t idx, const walberla::real_t& v) { temperature_[idx] = v; }
   
   const walberla::real_t& getHeatFlux(const size_t idx) const {return heatFlux_[idx];}
   walberla::real_t& getHeatFluxRef(const size_t idx) {return heatFlux_[idx];}
   void setHeatFlux(const size_t idx, const walberla::real_t& v) { heatFlux_[idx] = v; }
   

   /**
    * @brief creates a new particle and returns an iterator pointing to it
    *
    * \attention Use this function only if you know what you are doing!
    * Messing with the uid might break the simulation!
    * If you are unsure use create(bool) instead.
    * @param uid unique id of the particle to be created
    * @return iterator to the newly created particle
    */
   inline iterator create(const id_t& uid);
   inline iterator create(const bool global = false);
   inline iterator erase(iterator& it);
   /// Finds the entry corresponding to \p uid.
   /// \return iterator to the object or end iterator
   inline iterator find(const id_t& uid);
   inline void reserve(const size_t size);
   inline void clear();
   inline size_t size() const;

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
   std::vector<walberla::id_t> uid_ {};
   std::vector<walberla::mesa_pd::Vec3> position_ {};
   std::vector<walberla::real_t> interactionRadius_ {};
   std::vector<walberla::mesa_pd::data::particle_flags::FlagT> flags_ {};
   std::vector<int> owner_ {};
   std::vector<std::vector<int>> ghostOwners_ {};
   std::vector<size_t> shapeID_ {};
   std::vector<walberla::mesa_pd::Rot3> rotation_ {};
   std::vector<walberla::mesa_pd::Vec3> angularVelocity_ {};
   std::vector<walberla::mesa_pd::Vec3> torque_ {};
   std::vector<walberla::mesa_pd::Vec3> linearVelocity_ {};
   std::vector<walberla::real_t> invMass_ {};
   std::vector<walberla::mesa_pd::Vec3> force_ {};
   std::vector<walberla::mesa_pd::Vec3> oldForce_ {};
   std::vector<walberla::mesa_pd::Vec3> oldTorque_ {};
   std::vector<uint_t> type_ {};
   std::vector<int> nextParticle_ {};
   std::vector<std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>> oldContactHistory_ {};
   std::vector<std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>> newContactHistory_ {};
   std::vector<walberla::real_t> temperature_ {};
   std::vector<walberla::real_t> heatFlux_ {};
   std::unordered_map<id_t, size_t> uidToIdx_;
   static_assert(std::is_same<decltype(uid_)::value_type, id_t>::value,
                 "Property uid of type id_t is missing. This property is required!");
};
using Particle = ParticleStorage::Particle;

inline
ParticleStorage::Particle& ParticleStorage::Particle::operator=(const ParticleStorage::Particle& rhs)
{
   getUidRef() = rhs.getUid();
   getPositionRef() = rhs.getPosition();
   getInteractionRadiusRef() = rhs.getInteractionRadius();
   getFlagsRef() = rhs.getFlags();
   getOwnerRef() = rhs.getOwner();
   getGhostOwnersRef() = rhs.getGhostOwners();
   getShapeIDRef() = rhs.getShapeID();
   getRotationRef() = rhs.getRotation();
   getAngularVelocityRef() = rhs.getAngularVelocity();
   getTorqueRef() = rhs.getTorque();
   getLinearVelocityRef() = rhs.getLinearVelocity();
   getInvMassRef() = rhs.getInvMass();
   getForceRef() = rhs.getForce();
   getOldForceRef() = rhs.getOldForce();
   getOldTorqueRef() = rhs.getOldTorque();
   getTypeRef() = rhs.getType();
   getNextParticleRef() = rhs.getNextParticle();
   getOldContactHistoryRef() = rhs.getOldContactHistory();
   getNewContactHistoryRef() = rhs.getNewContactHistory();
   getTemperatureRef() = rhs.getTemperature();
   getHeatFluxRef() = rhs.getHeatFlux();
   return *this;
}

inline
ParticleStorage::Particle& ParticleStorage::Particle::operator=(ParticleStorage::Particle&& rhs)
{
   getUidRef() = std::move(rhs.getUidRef());
   getPositionRef() = std::move(rhs.getPositionRef());
   getInteractionRadiusRef() = std::move(rhs.getInteractionRadiusRef());
   getFlagsRef() = std::move(rhs.getFlagsRef());
   getOwnerRef() = std::move(rhs.getOwnerRef());
   getGhostOwnersRef() = std::move(rhs.getGhostOwnersRef());
   getShapeIDRef() = std::move(rhs.getShapeIDRef());
   getRotationRef() = std::move(rhs.getRotationRef());
   getAngularVelocityRef() = std::move(rhs.getAngularVelocityRef());
   getTorqueRef() = std::move(rhs.getTorqueRef());
   getLinearVelocityRef() = std::move(rhs.getLinearVelocityRef());
   getInvMassRef() = std::move(rhs.getInvMassRef());
   getForceRef() = std::move(rhs.getForceRef());
   getOldForceRef() = std::move(rhs.getOldForceRef());
   getOldTorqueRef() = std::move(rhs.getOldTorqueRef());
   getTypeRef() = std::move(rhs.getTypeRef());
   getNextParticleRef() = std::move(rhs.getNextParticleRef());
   getOldContactHistoryRef() = std::move(rhs.getOldContactHistoryRef());
   getNewContactHistoryRef() = std::move(rhs.getNewContactHistoryRef());
   getTemperatureRef() = std::move(rhs.getTemperatureRef());
   getHeatFluxRef() = std::move(rhs.getHeatFluxRef());
   return *this;
}

inline
std::ostream& operator<<( std::ostream& os, const ParticleStorage::Particle& p )
{
   os << "==========    ==========" << "\n" <<
         "idx                 : " << p.getIdx() << "\n" <<
         "uid                 : " << p.getUid() << "\n" <<
         "position            : " << p.getPosition() << "\n" <<
         "interactionRadius   : " << p.getInteractionRadius() << "\n" <<
         "flags               : " << p.getFlags() << "\n" <<
         "owner               : " << p.getOwner() << "\n" <<
         "ghostOwners         : " << p.getGhostOwners() << "\n" <<
         "shapeID             : " << p.getShapeID() << "\n" <<
         "rotation            : " << p.getRotation() << "\n" <<
         "angularVelocity     : " << p.getAngularVelocity() << "\n" <<
         "torque              : " << p.getTorque() << "\n" <<
         "linearVelocity      : " << p.getLinearVelocity() << "\n" <<
         "invMass             : " << p.getInvMass() << "\n" <<
         "force               : " << p.getForce() << "\n" <<
         "oldForce            : " << p.getOldForce() << "\n" <<
         "oldTorque           : " << p.getOldTorque() << "\n" <<
         "type                : " << p.getType() << "\n" <<
         "nextParticle        : " << p.getNextParticle() << "\n" <<
         "oldContactHistory   : " << p.getOldContactHistory() << "\n" <<
         "newContactHistory   : " << p.getNewContactHistory() << "\n" <<
         "temperature         : " << p.getTemperature() << "\n" <<
         "heatFlux            : " << p.getHeatFlux() << "\n" <<
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
   uid_.emplace_back(UniqueID<data::Particle>::invalidID());
   position_.emplace_back(real_t(0));
   interactionRadius_.emplace_back(real_t(0));
   flags_.emplace_back();
   owner_.emplace_back(-1);
   ghostOwners_.emplace_back();
   shapeID_.emplace_back();
   rotation_.emplace_back();
   angularVelocity_.emplace_back(real_t(0));
   torque_.emplace_back(real_t(0));
   linearVelocity_.emplace_back(real_t(0));
   invMass_.emplace_back(real_t(1));
   force_.emplace_back(real_t(0));
   oldForce_.emplace_back(real_t(0));
   oldTorque_.emplace_back(real_t(0));
   type_.emplace_back(0);
   nextParticle_.emplace_back(-1);
   oldContactHistory_.emplace_back();
   newContactHistory_.emplace_back();
   temperature_.emplace_back(real_t(0));
   heatFlux_.emplace_back(real_t(0));
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
   uid_.pop_back();
   position_.pop_back();
   interactionRadius_.pop_back();
   flags_.pop_back();
   owner_.pop_back();
   ghostOwners_.pop_back();
   shapeID_.pop_back();
   rotation_.pop_back();
   angularVelocity_.pop_back();
   torque_.pop_back();
   linearVelocity_.pop_back();
   invMass_.pop_back();
   force_.pop_back();
   oldForce_.pop_back();
   oldTorque_.pop_back();
   type_.pop_back();
   nextParticle_.pop_back();
   oldContactHistory_.pop_back();
   newContactHistory_.pop_back();
   temperature_.pop_back();
   heatFlux_.pop_back();
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
   WALBERLA_ASSERT_EQUAL(it->first, uid, "Lookup via uidToIdx map is not up to date!!!");
   return iterator(this, it->second);
}

inline void ParticleStorage::reserve(const size_t size)
{
   uid_.reserve(size);
   position_.reserve(size);
   interactionRadius_.reserve(size);
   flags_.reserve(size);
   owner_.reserve(size);
   ghostOwners_.reserve(size);
   shapeID_.reserve(size);
   rotation_.reserve(size);
   angularVelocity_.reserve(size);
   torque_.reserve(size);
   linearVelocity_.reserve(size);
   invMass_.reserve(size);
   force_.reserve(size);
   oldForce_.reserve(size);
   oldTorque_.reserve(size);
   type_.reserve(size);
   nextParticle_.reserve(size);
   oldContactHistory_.reserve(size);
   newContactHistory_.reserve(size);
   temperature_.reserve(size);
   heatFlux_.reserve(size);
}

inline void ParticleStorage::clear()
{
   uid_.clear();
   position_.clear();
   interactionRadius_.clear();
   flags_.clear();
   owner_.clear();
   ghostOwners_.clear();
   shapeID_.clear();
   rotation_.clear();
   angularVelocity_.clear();
   torque_.clear();
   linearVelocity_.clear();
   invMass_.clear();
   force_.clear();
   oldForce_.clear();
   oldTorque_.clear();
   type_.clear();
   nextParticle_.clear();
   oldContactHistory_.clear();
   newContactHistory_.clear();
   temperature_.clear();
   heatFlux_.clear();
   uidToIdx_.clear();
}

inline size_t ParticleStorage::size() const
{
   //WALBERLA_ASSERT_EQUAL( uid_.size(), uid.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), position.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), interactionRadius.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), flags.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), owner.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), ghostOwners.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), shapeID.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), rotation.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), angularVelocity.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), torque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), linearVelocity.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), invMass.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), force.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldForce.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldTorque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), type.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), nextParticle.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldContactHistory.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), newContactHistory.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), temperature.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), heatFlux.size() );
   return uid_.size();
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticle(const bool openmp,
                                             const Selector& selector,
                                             Accessor& acForPS,
                                             Func&& func, Args&&... args) 
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
   for (int64_t i = 0; i < int64_c(len); ++i)
   {
      if (selector(uint64_c(i), acForPS))
         func( uint64_c(i), std::forward<Args>(args)... );
   }
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticle(const bool openmp,
                                             const Selector& selector,
                                             Accessor& acForPS,
                                             Func&& func, Args&&... args) const
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
   for (int64_t i = 0; i < int64_c(len); ++i)
   {
      if (selector(uint64_c(i), acForPS))
         func( uint64_c(i), std::forward<Args>(args)... );
   }
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePair(const bool openmp,
                                                 const Selector& selector,
                                                 Accessor& acForPS,
                                                 Func&& func, Args&&... args) 
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
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
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePair(const bool openmp,
                                                 const Selector& selector,
                                                 Accessor& acForPS,
                                                 Func&& func, Args&&... args) const
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
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
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePairHalf(const bool openmp,
                                                     const Selector& selector,
                                                     Accessor& acForPS,
                                                     Func&& func, Args&&... args) 
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
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
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticlePairHalf(const bool openmp,
                                                     const Selector& selector,
                                                     Accessor& acForPS,
                                                     Func&& func, Args&&... args) const
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if (openmp)
   #endif
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
///Predicate that selects a certain property from a Particle
class SelectParticleUid
{
public:
   using return_type = walberla::id_t;
   walberla::id_t& operator()(data::Particle& p) const {return p.getUidRef();}
   walberla::id_t& operator()(data::Particle&& p) const {return p.getUidRef();}
   const walberla::id_t& operator()(const data::Particle& p) const {return p.getUid();}
};
///Predicate that selects a certain property from a Particle
class SelectParticlePosition
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getPositionRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getPositionRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getPosition();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleInteractionRadius
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getInteractionRadiusRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getInteractionRadiusRef();}
   const walberla::real_t& operator()(const data::Particle& p) const {return p.getInteractionRadius();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleFlags
{
public:
   using return_type = walberla::mesa_pd::data::particle_flags::FlagT;
   walberla::mesa_pd::data::particle_flags::FlagT& operator()(data::Particle& p) const {return p.getFlagsRef();}
   walberla::mesa_pd::data::particle_flags::FlagT& operator()(data::Particle&& p) const {return p.getFlagsRef();}
   const walberla::mesa_pd::data::particle_flags::FlagT& operator()(const data::Particle& p) const {return p.getFlags();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOwner
{
public:
   using return_type = int;
   int& operator()(data::Particle& p) const {return p.getOwnerRef();}
   int& operator()(data::Particle&& p) const {return p.getOwnerRef();}
   const int& operator()(const data::Particle& p) const {return p.getOwner();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleGhostOwners
{
public:
   using return_type = std::vector<int>;
   std::vector<int>& operator()(data::Particle& p) const {return p.getGhostOwnersRef();}
   std::vector<int>& operator()(data::Particle&& p) const {return p.getGhostOwnersRef();}
   const std::vector<int>& operator()(const data::Particle& p) const {return p.getGhostOwners();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleShapeID
{
public:
   using return_type = size_t;
   size_t& operator()(data::Particle& p) const {return p.getShapeIDRef();}
   size_t& operator()(data::Particle&& p) const {return p.getShapeIDRef();}
   const size_t& operator()(const data::Particle& p) const {return p.getShapeID();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleRotation
{
public:
   using return_type = walberla::mesa_pd::Rot3;
   walberla::mesa_pd::Rot3& operator()(data::Particle& p) const {return p.getRotationRef();}
   walberla::mesa_pd::Rot3& operator()(data::Particle&& p) const {return p.getRotationRef();}
   const walberla::mesa_pd::Rot3& operator()(const data::Particle& p) const {return p.getRotation();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleAngularVelocity
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getAngularVelocityRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getAngularVelocityRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getAngularVelocity();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getTorqueRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleLinearVelocity
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getLinearVelocityRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getLinearVelocityRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getLinearVelocity();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleInvMass
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getInvMassRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getInvMassRef();}
   const walberla::real_t& operator()(const data::Particle& p) const {return p.getInvMass();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getForceRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldForceRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getOldForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldTorqueRef();}
   const walberla::mesa_pd::Vec3& operator()(const data::Particle& p) const {return p.getOldTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleType
{
public:
   using return_type = uint_t;
   uint_t& operator()(data::Particle& p) const {return p.getTypeRef();}
   uint_t& operator()(data::Particle&& p) const {return p.getTypeRef();}
   const uint_t& operator()(const data::Particle& p) const {return p.getType();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleNextParticle
{
public:
   using return_type = int;
   int& operator()(data::Particle& p) const {return p.getNextParticleRef();}
   int& operator()(data::Particle&& p) const {return p.getNextParticleRef();}
   const int& operator()(const data::Particle& p) const {return p.getNextParticle();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldContactHistory
{
public:
   using return_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle& p) const {return p.getOldContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle&& p) const {return p.getOldContactHistoryRef();}
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(const data::Particle& p) const {return p.getOldContactHistory();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleNewContactHistory
{
public:
   using return_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle& p) const {return p.getNewContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle&& p) const {return p.getNewContactHistoryRef();}
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(const data::Particle& p) const {return p.getNewContactHistory();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleTemperature
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getTemperatureRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getTemperatureRef();}
   const walberla::real_t& operator()(const data::Particle& p) const {return p.getTemperature();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleHeatFlux
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getHeatFluxRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getHeatFluxRef();}
   const walberla::real_t& operator()(const data::Particle& p) const {return p.getHeatFlux();}
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla