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
#include <mesa_pd/data/Flags.h>
#include <blockforest/BlockForest.h>
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

      
      using uid_type = walberla::id_t;
      using position_type = walberla::mesa_pd::Vec3;
      using interactionRadius_type = walberla::real_t;
      using flags_type = walberla::mesa_pd::data::particle_flags::FlagT;
      using owner_type = int;
      using ghostOwners_type = std::unordered_set<walberla::mpi::MPIRank>;
      using linearVelocity_type = walberla::mesa_pd::Vec3;
      using invMass_type = walberla::real_t;
      using force_type = walberla::mesa_pd::Vec3;
      using oldForce_type = walberla::mesa_pd::Vec3;
      using shapeID_type = size_t;
      using rotation_type = walberla::mesa_pd::Rot3;
      using angularVelocity_type = walberla::mesa_pd::Vec3;
      using torque_type = walberla::mesa_pd::Vec3;
      using oldTorque_type = walberla::mesa_pd::Vec3;
      using radiusAtTemperature_type = walberla::real_t;
      using currentBlock_type = blockforest::BlockID;
      using type_type = uint_t;
      using nextParticle_type = int;
      using oldContactHistory_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
      using newContactHistory_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
      using temperature_type = walberla::real_t;
      using heatFlux_type = walberla::real_t;
      using dv_type = walberla::mesa_pd::Vec3;
      using dw_type = walberla::mesa_pd::Vec3;
      using hydrodynamicForce_type = walberla::mesa_pd::Vec3;
      using hydrodynamicTorque_type = walberla::mesa_pd::Vec3;
      using oldHydrodynamicForce_type = walberla::mesa_pd::Vec3;
      using oldHydrodynamicTorque_type = walberla::mesa_pd::Vec3;
      using neighborState_type = std::unordered_set<walberla::mpi::MPIRank>;

      
      uid_type const & getUid() const {return storage_.getUid(i_);}
      uid_type& getUidRef() {return storage_.getUidRef(i_);}
      void setUid(uid_type const & v) { storage_.setUid(i_, v);}
      
      position_type const & getPosition() const {return storage_.getPosition(i_);}
      position_type& getPositionRef() {return storage_.getPositionRef(i_);}
      void setPosition(position_type const & v) { storage_.setPosition(i_, v);}
      
      interactionRadius_type const & getInteractionRadius() const {return storage_.getInteractionRadius(i_);}
      interactionRadius_type& getInteractionRadiusRef() {return storage_.getInteractionRadiusRef(i_);}
      void setInteractionRadius(interactionRadius_type const & v) { storage_.setInteractionRadius(i_, v);}
      
      flags_type const & getFlags() const {return storage_.getFlags(i_);}
      flags_type& getFlagsRef() {return storage_.getFlagsRef(i_);}
      void setFlags(flags_type const & v) { storage_.setFlags(i_, v);}
      
      owner_type const & getOwner() const {return storage_.getOwner(i_);}
      owner_type& getOwnerRef() {return storage_.getOwnerRef(i_);}
      void setOwner(owner_type const & v) { storage_.setOwner(i_, v);}
      
      ghostOwners_type const & getGhostOwners() const {return storage_.getGhostOwners(i_);}
      ghostOwners_type& getGhostOwnersRef() {return storage_.getGhostOwnersRef(i_);}
      void setGhostOwners(ghostOwners_type const & v) { storage_.setGhostOwners(i_, v);}
      
      linearVelocity_type const & getLinearVelocity() const {return storage_.getLinearVelocity(i_);}
      linearVelocity_type& getLinearVelocityRef() {return storage_.getLinearVelocityRef(i_);}
      void setLinearVelocity(linearVelocity_type const & v) { storage_.setLinearVelocity(i_, v);}
      
      invMass_type const & getInvMass() const {return storage_.getInvMass(i_);}
      invMass_type& getInvMassRef() {return storage_.getInvMassRef(i_);}
      void setInvMass(invMass_type const & v) { storage_.setInvMass(i_, v);}
      
      force_type const & getForce() const {return storage_.getForce(i_);}
      force_type& getForceRef() {return storage_.getForceRef(i_);}
      void setForce(force_type const & v) { storage_.setForce(i_, v);}
      
      oldForce_type const & getOldForce() const {return storage_.getOldForce(i_);}
      oldForce_type& getOldForceRef() {return storage_.getOldForceRef(i_);}
      void setOldForce(oldForce_type const & v) { storage_.setOldForce(i_, v);}
      
      shapeID_type const & getShapeID() const {return storage_.getShapeID(i_);}
      shapeID_type& getShapeIDRef() {return storage_.getShapeIDRef(i_);}
      void setShapeID(shapeID_type const & v) { storage_.setShapeID(i_, v);}
      
      rotation_type const & getRotation() const {return storage_.getRotation(i_);}
      rotation_type& getRotationRef() {return storage_.getRotationRef(i_);}
      void setRotation(rotation_type const & v) { storage_.setRotation(i_, v);}
      
      angularVelocity_type const & getAngularVelocity() const {return storage_.getAngularVelocity(i_);}
      angularVelocity_type& getAngularVelocityRef() {return storage_.getAngularVelocityRef(i_);}
      void setAngularVelocity(angularVelocity_type const & v) { storage_.setAngularVelocity(i_, v);}
      
      torque_type const & getTorque() const {return storage_.getTorque(i_);}
      torque_type& getTorqueRef() {return storage_.getTorqueRef(i_);}
      void setTorque(torque_type const & v) { storage_.setTorque(i_, v);}
      
      oldTorque_type const & getOldTorque() const {return storage_.getOldTorque(i_);}
      oldTorque_type& getOldTorqueRef() {return storage_.getOldTorqueRef(i_);}
      void setOldTorque(oldTorque_type const & v) { storage_.setOldTorque(i_, v);}
      
      radiusAtTemperature_type const & getRadiusAtTemperature() const {return storage_.getRadiusAtTemperature(i_);}
      radiusAtTemperature_type& getRadiusAtTemperatureRef() {return storage_.getRadiusAtTemperatureRef(i_);}
      void setRadiusAtTemperature(radiusAtTemperature_type const & v) { storage_.setRadiusAtTemperature(i_, v);}
      
      currentBlock_type const & getCurrentBlock() const {return storage_.getCurrentBlock(i_);}
      currentBlock_type& getCurrentBlockRef() {return storage_.getCurrentBlockRef(i_);}
      void setCurrentBlock(currentBlock_type const & v) { storage_.setCurrentBlock(i_, v);}
      
      type_type const & getType() const {return storage_.getType(i_);}
      type_type& getTypeRef() {return storage_.getTypeRef(i_);}
      void setType(type_type const & v) { storage_.setType(i_, v);}
      
      nextParticle_type const & getNextParticle() const {return storage_.getNextParticle(i_);}
      nextParticle_type& getNextParticleRef() {return storage_.getNextParticleRef(i_);}
      void setNextParticle(nextParticle_type const & v) { storage_.setNextParticle(i_, v);}
      
      oldContactHistory_type const & getOldContactHistory() const {return storage_.getOldContactHistory(i_);}
      oldContactHistory_type& getOldContactHistoryRef() {return storage_.getOldContactHistoryRef(i_);}
      void setOldContactHistory(oldContactHistory_type const & v) { storage_.setOldContactHistory(i_, v);}
      
      newContactHistory_type const & getNewContactHistory() const {return storage_.getNewContactHistory(i_);}
      newContactHistory_type& getNewContactHistoryRef() {return storage_.getNewContactHistoryRef(i_);}
      void setNewContactHistory(newContactHistory_type const & v) { storage_.setNewContactHistory(i_, v);}
      
      temperature_type const & getTemperature() const {return storage_.getTemperature(i_);}
      temperature_type& getTemperatureRef() {return storage_.getTemperatureRef(i_);}
      void setTemperature(temperature_type const & v) { storage_.setTemperature(i_, v);}
      
      heatFlux_type const & getHeatFlux() const {return storage_.getHeatFlux(i_);}
      heatFlux_type& getHeatFluxRef() {return storage_.getHeatFluxRef(i_);}
      void setHeatFlux(heatFlux_type const & v) { storage_.setHeatFlux(i_, v);}
      
      dv_type const & getDv() const {return storage_.getDv(i_);}
      dv_type& getDvRef() {return storage_.getDvRef(i_);}
      void setDv(dv_type const & v) { storage_.setDv(i_, v);}
      
      dw_type const & getDw() const {return storage_.getDw(i_);}
      dw_type& getDwRef() {return storage_.getDwRef(i_);}
      void setDw(dw_type const & v) { storage_.setDw(i_, v);}
      
      hydrodynamicForce_type const & getHydrodynamicForce() const {return storage_.getHydrodynamicForce(i_);}
      hydrodynamicForce_type& getHydrodynamicForceRef() {return storage_.getHydrodynamicForceRef(i_);}
      void setHydrodynamicForce(hydrodynamicForce_type const & v) { storage_.setHydrodynamicForce(i_, v);}
      
      hydrodynamicTorque_type const & getHydrodynamicTorque() const {return storage_.getHydrodynamicTorque(i_);}
      hydrodynamicTorque_type& getHydrodynamicTorqueRef() {return storage_.getHydrodynamicTorqueRef(i_);}
      void setHydrodynamicTorque(hydrodynamicTorque_type const & v) { storage_.setHydrodynamicTorque(i_, v);}
      
      oldHydrodynamicForce_type const & getOldHydrodynamicForce() const {return storage_.getOldHydrodynamicForce(i_);}
      oldHydrodynamicForce_type& getOldHydrodynamicForceRef() {return storage_.getOldHydrodynamicForceRef(i_);}
      void setOldHydrodynamicForce(oldHydrodynamicForce_type const & v) { storage_.setOldHydrodynamicForce(i_, v);}
      
      oldHydrodynamicTorque_type const & getOldHydrodynamicTorque() const {return storage_.getOldHydrodynamicTorque(i_);}
      oldHydrodynamicTorque_type& getOldHydrodynamicTorqueRef() {return storage_.getOldHydrodynamicTorqueRef(i_);}
      void setOldHydrodynamicTorque(oldHydrodynamicTorque_type const & v) { storage_.setOldHydrodynamicTorque(i_, v);}
      
      neighborState_type const & getNeighborState() const {return storage_.getNeighborState(i_);}
      neighborState_type& getNeighborStateRef() {return storage_.getNeighborStateRef(i_);}
      void setNeighborState(neighborState_type const & v) { storage_.setNeighborState(i_, v);}
      

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

   
   using uid_type = walberla::id_t;
   using position_type = walberla::mesa_pd::Vec3;
   using interactionRadius_type = walberla::real_t;
   using flags_type = walberla::mesa_pd::data::particle_flags::FlagT;
   using owner_type = int;
   using ghostOwners_type = std::unordered_set<walberla::mpi::MPIRank>;
   using linearVelocity_type = walberla::mesa_pd::Vec3;
   using invMass_type = walberla::real_t;
   using force_type = walberla::mesa_pd::Vec3;
   using oldForce_type = walberla::mesa_pd::Vec3;
   using shapeID_type = size_t;
   using rotation_type = walberla::mesa_pd::Rot3;
   using angularVelocity_type = walberla::mesa_pd::Vec3;
   using torque_type = walberla::mesa_pd::Vec3;
   using oldTorque_type = walberla::mesa_pd::Vec3;
   using radiusAtTemperature_type = walberla::real_t;
   using currentBlock_type = blockforest::BlockID;
   using type_type = uint_t;
   using nextParticle_type = int;
   using oldContactHistory_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   using newContactHistory_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   using temperature_type = walberla::real_t;
   using heatFlux_type = walberla::real_t;
   using dv_type = walberla::mesa_pd::Vec3;
   using dw_type = walberla::mesa_pd::Vec3;
   using hydrodynamicForce_type = walberla::mesa_pd::Vec3;
   using hydrodynamicTorque_type = walberla::mesa_pd::Vec3;
   using oldHydrodynamicForce_type = walberla::mesa_pd::Vec3;
   using oldHydrodynamicTorque_type = walberla::mesa_pd::Vec3;
   using neighborState_type = std::unordered_set<walberla::mpi::MPIRank>;

   
   uid_type const & getUid(const size_t idx) const {return uid_[idx];}
   uid_type& getUidRef(const size_t idx) {return uid_[idx];}
   void setUid(const size_t idx, uid_type const & v) { uid_[idx] = v; }
   
   position_type const & getPosition(const size_t idx) const {return position_[idx];}
   position_type& getPositionRef(const size_t idx) {return position_[idx];}
   void setPosition(const size_t idx, position_type const & v) { position_[idx] = v; }
   
   interactionRadius_type const & getInteractionRadius(const size_t idx) const {return interactionRadius_[idx];}
   interactionRadius_type& getInteractionRadiusRef(const size_t idx) {return interactionRadius_[idx];}
   void setInteractionRadius(const size_t idx, interactionRadius_type const & v) { interactionRadius_[idx] = v; }
   
   flags_type const & getFlags(const size_t idx) const {return flags_[idx];}
   flags_type& getFlagsRef(const size_t idx) {return flags_[idx];}
   void setFlags(const size_t idx, flags_type const & v) { flags_[idx] = v; }
   
   owner_type const & getOwner(const size_t idx) const {return owner_[idx];}
   owner_type& getOwnerRef(const size_t idx) {return owner_[idx];}
   void setOwner(const size_t idx, owner_type const & v) { owner_[idx] = v; }
   
   ghostOwners_type const & getGhostOwners(const size_t idx) const {return ghostOwners_[idx];}
   ghostOwners_type& getGhostOwnersRef(const size_t idx) {return ghostOwners_[idx];}
   void setGhostOwners(const size_t idx, ghostOwners_type const & v) { ghostOwners_[idx] = v; }
   
   linearVelocity_type const & getLinearVelocity(const size_t idx) const {return linearVelocity_[idx];}
   linearVelocity_type& getLinearVelocityRef(const size_t idx) {return linearVelocity_[idx];}
   void setLinearVelocity(const size_t idx, linearVelocity_type const & v) { linearVelocity_[idx] = v; }
   
   invMass_type const & getInvMass(const size_t idx) const {return invMass_[idx];}
   invMass_type& getInvMassRef(const size_t idx) {return invMass_[idx];}
   void setInvMass(const size_t idx, invMass_type const & v) { invMass_[idx] = v; }
   
   force_type const & getForce(const size_t idx) const {return force_[idx];}
   force_type& getForceRef(const size_t idx) {return force_[idx];}
   void setForce(const size_t idx, force_type const & v) { force_[idx] = v; }
   
   oldForce_type const & getOldForce(const size_t idx) const {return oldForce_[idx];}
   oldForce_type& getOldForceRef(const size_t idx) {return oldForce_[idx];}
   void setOldForce(const size_t idx, oldForce_type const & v) { oldForce_[idx] = v; }
   
   shapeID_type const & getShapeID(const size_t idx) const {return shapeID_[idx];}
   shapeID_type& getShapeIDRef(const size_t idx) {return shapeID_[idx];}
   void setShapeID(const size_t idx, shapeID_type const & v) { shapeID_[idx] = v; }
   
   rotation_type const & getRotation(const size_t idx) const {return rotation_[idx];}
   rotation_type& getRotationRef(const size_t idx) {return rotation_[idx];}
   void setRotation(const size_t idx, rotation_type const & v) { rotation_[idx] = v; }
   
   angularVelocity_type const & getAngularVelocity(const size_t idx) const {return angularVelocity_[idx];}
   angularVelocity_type& getAngularVelocityRef(const size_t idx) {return angularVelocity_[idx];}
   void setAngularVelocity(const size_t idx, angularVelocity_type const & v) { angularVelocity_[idx] = v; }
   
   torque_type const & getTorque(const size_t idx) const {return torque_[idx];}
   torque_type& getTorqueRef(const size_t idx) {return torque_[idx];}
   void setTorque(const size_t idx, torque_type const & v) { torque_[idx] = v; }
   
   oldTorque_type const & getOldTorque(const size_t idx) const {return oldTorque_[idx];}
   oldTorque_type& getOldTorqueRef(const size_t idx) {return oldTorque_[idx];}
   void setOldTorque(const size_t idx, oldTorque_type const & v) { oldTorque_[idx] = v; }
   
   radiusAtTemperature_type const & getRadiusAtTemperature(const size_t idx) const {return radiusAtTemperature_[idx];}
   radiusAtTemperature_type& getRadiusAtTemperatureRef(const size_t idx) {return radiusAtTemperature_[idx];}
   void setRadiusAtTemperature(const size_t idx, radiusAtTemperature_type const & v) { radiusAtTemperature_[idx] = v; }
   
   currentBlock_type const & getCurrentBlock(const size_t idx) const {return currentBlock_[idx];}
   currentBlock_type& getCurrentBlockRef(const size_t idx) {return currentBlock_[idx];}
   void setCurrentBlock(const size_t idx, currentBlock_type const & v) { currentBlock_[idx] = v; }
   
   type_type const & getType(const size_t idx) const {return type_[idx];}
   type_type& getTypeRef(const size_t idx) {return type_[idx];}
   void setType(const size_t idx, type_type const & v) { type_[idx] = v; }
   
   nextParticle_type const & getNextParticle(const size_t idx) const {return nextParticle_[idx];}
   nextParticle_type& getNextParticleRef(const size_t idx) {return nextParticle_[idx];}
   void setNextParticle(const size_t idx, nextParticle_type const & v) { nextParticle_[idx] = v; }
   
   oldContactHistory_type const & getOldContactHistory(const size_t idx) const {return oldContactHistory_[idx];}
   oldContactHistory_type& getOldContactHistoryRef(const size_t idx) {return oldContactHistory_[idx];}
   void setOldContactHistory(const size_t idx, oldContactHistory_type const & v) { oldContactHistory_[idx] = v; }
   
   newContactHistory_type const & getNewContactHistory(const size_t idx) const {return newContactHistory_[idx];}
   newContactHistory_type& getNewContactHistoryRef(const size_t idx) {return newContactHistory_[idx];}
   void setNewContactHistory(const size_t idx, newContactHistory_type const & v) { newContactHistory_[idx] = v; }
   
   temperature_type const & getTemperature(const size_t idx) const {return temperature_[idx];}
   temperature_type& getTemperatureRef(const size_t idx) {return temperature_[idx];}
   void setTemperature(const size_t idx, temperature_type const & v) { temperature_[idx] = v; }
   
   heatFlux_type const & getHeatFlux(const size_t idx) const {return heatFlux_[idx];}
   heatFlux_type& getHeatFluxRef(const size_t idx) {return heatFlux_[idx];}
   void setHeatFlux(const size_t idx, heatFlux_type const & v) { heatFlux_[idx] = v; }
   
   dv_type const & getDv(const size_t idx) const {return dv_[idx];}
   dv_type& getDvRef(const size_t idx) {return dv_[idx];}
   void setDv(const size_t idx, dv_type const & v) { dv_[idx] = v; }
   
   dw_type const & getDw(const size_t idx) const {return dw_[idx];}
   dw_type& getDwRef(const size_t idx) {return dw_[idx];}
   void setDw(const size_t idx, dw_type const & v) { dw_[idx] = v; }
   
   hydrodynamicForce_type const & getHydrodynamicForce(const size_t idx) const {return hydrodynamicForce_[idx];}
   hydrodynamicForce_type& getHydrodynamicForceRef(const size_t idx) {return hydrodynamicForce_[idx];}
   void setHydrodynamicForce(const size_t idx, hydrodynamicForce_type const & v) { hydrodynamicForce_[idx] = v; }
   
   hydrodynamicTorque_type const & getHydrodynamicTorque(const size_t idx) const {return hydrodynamicTorque_[idx];}
   hydrodynamicTorque_type& getHydrodynamicTorqueRef(const size_t idx) {return hydrodynamicTorque_[idx];}
   void setHydrodynamicTorque(const size_t idx, hydrodynamicTorque_type const & v) { hydrodynamicTorque_[idx] = v; }
   
   oldHydrodynamicForce_type const & getOldHydrodynamicForce(const size_t idx) const {return oldHydrodynamicForce_[idx];}
   oldHydrodynamicForce_type& getOldHydrodynamicForceRef(const size_t idx) {return oldHydrodynamicForce_[idx];}
   void setOldHydrodynamicForce(const size_t idx, oldHydrodynamicForce_type const & v) { oldHydrodynamicForce_[idx] = v; }
   
   oldHydrodynamicTorque_type const & getOldHydrodynamicTorque(const size_t idx) const {return oldHydrodynamicTorque_[idx];}
   oldHydrodynamicTorque_type& getOldHydrodynamicTorqueRef(const size_t idx) {return oldHydrodynamicTorque_[idx];}
   void setOldHydrodynamicTorque(const size_t idx, oldHydrodynamicTorque_type const & v) { oldHydrodynamicTorque_[idx] = v; }
   
   neighborState_type const & getNeighborState(const size_t idx) const {return neighborState_[idx];}
   neighborState_type& getNeighborStateRef(const size_t idx) {return neighborState_[idx];}
   void setNeighborState(const size_t idx, neighborState_type const & v) { neighborState_[idx] = v; }
   

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
   std::vector<uid_type> uid_ {};
   std::vector<position_type> position_ {};
   std::vector<interactionRadius_type> interactionRadius_ {};
   std::vector<flags_type> flags_ {};
   std::vector<owner_type> owner_ {};
   std::vector<ghostOwners_type> ghostOwners_ {};
   std::vector<linearVelocity_type> linearVelocity_ {};
   std::vector<invMass_type> invMass_ {};
   std::vector<force_type> force_ {};
   std::vector<oldForce_type> oldForce_ {};
   std::vector<shapeID_type> shapeID_ {};
   std::vector<rotation_type> rotation_ {};
   std::vector<angularVelocity_type> angularVelocity_ {};
   std::vector<torque_type> torque_ {};
   std::vector<oldTorque_type> oldTorque_ {};
   std::vector<radiusAtTemperature_type> radiusAtTemperature_ {};
   std::vector<currentBlock_type> currentBlock_ {};
   std::vector<type_type> type_ {};
   std::vector<nextParticle_type> nextParticle_ {};
   std::vector<oldContactHistory_type> oldContactHistory_ {};
   std::vector<newContactHistory_type> newContactHistory_ {};
   std::vector<temperature_type> temperature_ {};
   std::vector<heatFlux_type> heatFlux_ {};
   std::vector<dv_type> dv_ {};
   std::vector<dw_type> dw_ {};
   std::vector<hydrodynamicForce_type> hydrodynamicForce_ {};
   std::vector<hydrodynamicTorque_type> hydrodynamicTorque_ {};
   std::vector<oldHydrodynamicForce_type> oldHydrodynamicForce_ {};
   std::vector<oldHydrodynamicTorque_type> oldHydrodynamicTorque_ {};
   std::vector<neighborState_type> neighborState_ {};
   std::unordered_map<uid_type, size_t> uidToIdx_;
   static_assert(std::is_same<uid_type, id_t>::value,
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
   getLinearVelocityRef() = rhs.getLinearVelocity();
   getInvMassRef() = rhs.getInvMass();
   getForceRef() = rhs.getForce();
   getOldForceRef() = rhs.getOldForce();
   getShapeIDRef() = rhs.getShapeID();
   getRotationRef() = rhs.getRotation();
   getAngularVelocityRef() = rhs.getAngularVelocity();
   getTorqueRef() = rhs.getTorque();
   getOldTorqueRef() = rhs.getOldTorque();
   getRadiusAtTemperatureRef() = rhs.getRadiusAtTemperature();
   getCurrentBlockRef() = rhs.getCurrentBlock();
   getTypeRef() = rhs.getType();
   getNextParticleRef() = rhs.getNextParticle();
   getOldContactHistoryRef() = rhs.getOldContactHistory();
   getNewContactHistoryRef() = rhs.getNewContactHistory();
   getTemperatureRef() = rhs.getTemperature();
   getHeatFluxRef() = rhs.getHeatFlux();
   getDvRef() = rhs.getDv();
   getDwRef() = rhs.getDw();
   getHydrodynamicForceRef() = rhs.getHydrodynamicForce();
   getHydrodynamicTorqueRef() = rhs.getHydrodynamicTorque();
   getOldHydrodynamicForceRef() = rhs.getOldHydrodynamicForce();
   getOldHydrodynamicTorqueRef() = rhs.getOldHydrodynamicTorque();
   getNeighborStateRef() = rhs.getNeighborState();
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
   getLinearVelocityRef() = std::move(rhs.getLinearVelocityRef());
   getInvMassRef() = std::move(rhs.getInvMassRef());
   getForceRef() = std::move(rhs.getForceRef());
   getOldForceRef() = std::move(rhs.getOldForceRef());
   getShapeIDRef() = std::move(rhs.getShapeIDRef());
   getRotationRef() = std::move(rhs.getRotationRef());
   getAngularVelocityRef() = std::move(rhs.getAngularVelocityRef());
   getTorqueRef() = std::move(rhs.getTorqueRef());
   getOldTorqueRef() = std::move(rhs.getOldTorqueRef());
   getRadiusAtTemperatureRef() = std::move(rhs.getRadiusAtTemperatureRef());
   getCurrentBlockRef() = std::move(rhs.getCurrentBlockRef());
   getTypeRef() = std::move(rhs.getTypeRef());
   getNextParticleRef() = std::move(rhs.getNextParticleRef());
   getOldContactHistoryRef() = std::move(rhs.getOldContactHistoryRef());
   getNewContactHistoryRef() = std::move(rhs.getNewContactHistoryRef());
   getTemperatureRef() = std::move(rhs.getTemperatureRef());
   getHeatFluxRef() = std::move(rhs.getHeatFluxRef());
   getDvRef() = std::move(rhs.getDvRef());
   getDwRef() = std::move(rhs.getDwRef());
   getHydrodynamicForceRef() = std::move(rhs.getHydrodynamicForceRef());
   getHydrodynamicTorqueRef() = std::move(rhs.getHydrodynamicTorqueRef());
   getOldHydrodynamicForceRef() = std::move(rhs.getOldHydrodynamicForceRef());
   getOldHydrodynamicTorqueRef() = std::move(rhs.getOldHydrodynamicTorqueRef());
   getNeighborStateRef() = std::move(rhs.getNeighborStateRef());
   return *this;
}

inline
void swap(ParticleStorage::Particle lhs, ParticleStorage::Particle rhs)
{
   if (lhs.i_ == rhs.i_) return;
   std::swap(lhs.getUidRef(), rhs.getUidRef());
   std::swap(lhs.getPositionRef(), rhs.getPositionRef());
   std::swap(lhs.getInteractionRadiusRef(), rhs.getInteractionRadiusRef());
   std::swap(lhs.getFlagsRef(), rhs.getFlagsRef());
   std::swap(lhs.getOwnerRef(), rhs.getOwnerRef());
   std::swap(lhs.getGhostOwnersRef(), rhs.getGhostOwnersRef());
   std::swap(lhs.getLinearVelocityRef(), rhs.getLinearVelocityRef());
   std::swap(lhs.getInvMassRef(), rhs.getInvMassRef());
   std::swap(lhs.getForceRef(), rhs.getForceRef());
   std::swap(lhs.getOldForceRef(), rhs.getOldForceRef());
   std::swap(lhs.getShapeIDRef(), rhs.getShapeIDRef());
   std::swap(lhs.getRotationRef(), rhs.getRotationRef());
   std::swap(lhs.getAngularVelocityRef(), rhs.getAngularVelocityRef());
   std::swap(lhs.getTorqueRef(), rhs.getTorqueRef());
   std::swap(lhs.getOldTorqueRef(), rhs.getOldTorqueRef());
   std::swap(lhs.getRadiusAtTemperatureRef(), rhs.getRadiusAtTemperatureRef());
   std::swap(lhs.getCurrentBlockRef(), rhs.getCurrentBlockRef());
   std::swap(lhs.getTypeRef(), rhs.getTypeRef());
   std::swap(lhs.getNextParticleRef(), rhs.getNextParticleRef());
   std::swap(lhs.getOldContactHistoryRef(), rhs.getOldContactHistoryRef());
   std::swap(lhs.getNewContactHistoryRef(), rhs.getNewContactHistoryRef());
   std::swap(lhs.getTemperatureRef(), rhs.getTemperatureRef());
   std::swap(lhs.getHeatFluxRef(), rhs.getHeatFluxRef());
   std::swap(lhs.getDvRef(), rhs.getDvRef());
   std::swap(lhs.getDwRef(), rhs.getDwRef());
   std::swap(lhs.getHydrodynamicForceRef(), rhs.getHydrodynamicForceRef());
   std::swap(lhs.getHydrodynamicTorqueRef(), rhs.getHydrodynamicTorqueRef());
   std::swap(lhs.getOldHydrodynamicForceRef(), rhs.getOldHydrodynamicForceRef());
   std::swap(lhs.getOldHydrodynamicTorqueRef(), rhs.getOldHydrodynamicTorqueRef());
   std::swap(lhs.getNeighborStateRef(), rhs.getNeighborStateRef());
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
         "linearVelocity      : " << p.getLinearVelocity() << "\n" <<
         "invMass             : " << p.getInvMass() << "\n" <<
         "force               : " << p.getForce() << "\n" <<
         "oldForce            : " << p.getOldForce() << "\n" <<
         "shapeID             : " << p.getShapeID() << "\n" <<
         "rotation            : " << p.getRotation() << "\n" <<
         "angularVelocity     : " << p.getAngularVelocity() << "\n" <<
         "torque              : " << p.getTorque() << "\n" <<
         "oldTorque           : " << p.getOldTorque() << "\n" <<
         "radiusAtTemperature : " << p.getRadiusAtTemperature() << "\n" <<
         "currentBlock        : " << p.getCurrentBlock() << "\n" <<
         "type                : " << p.getType() << "\n" <<
         "nextParticle        : " << p.getNextParticle() << "\n" <<
         "oldContactHistory   : " << p.getOldContactHistory() << "\n" <<
         "newContactHistory   : " << p.getNewContactHistory() << "\n" <<
         "temperature         : " << p.getTemperature() << "\n" <<
         "heatFlux            : " << p.getHeatFlux() << "\n" <<
         "dv                  : " << p.getDv() << "\n" <<
         "dw                  : " << p.getDw() << "\n" <<
         "hydrodynamicForce   : " << p.getHydrodynamicForce() << "\n" <<
         "hydrodynamicTorque  : " << p.getHydrodynamicTorque() << "\n" <<
         "oldHydrodynamicForce: " << p.getOldHydrodynamicForce() << "\n" <<
         "oldHydrodynamicTorque: " << p.getOldHydrodynamicTorque() << "\n" <<
         "neighborState       : " << p.getNeighborState() << "\n" <<
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
   linearVelocity_.emplace_back(real_t(0));
   invMass_.emplace_back(real_t(1));
   force_.emplace_back(real_t(0));
   oldForce_.emplace_back(real_t(0));
   shapeID_.emplace_back();
   rotation_.emplace_back();
   angularVelocity_.emplace_back(real_t(0));
   torque_.emplace_back(real_t(0));
   oldTorque_.emplace_back(real_t(0));
   radiusAtTemperature_.emplace_back(real_t(0));
   currentBlock_.emplace_back();
   type_.emplace_back(0);
   nextParticle_.emplace_back(-1);
   oldContactHistory_.emplace_back();
   newContactHistory_.emplace_back();
   temperature_.emplace_back(real_t(0));
   heatFlux_.emplace_back(real_t(0));
   dv_.emplace_back(real_t(0));
   dw_.emplace_back(real_t(0));
   hydrodynamicForce_.emplace_back(real_t(0));
   hydrodynamicTorque_.emplace_back(real_t(0));
   oldHydrodynamicForce_.emplace_back(real_t(0));
   oldHydrodynamicTorque_.emplace_back(real_t(0));
   neighborState_.emplace_back();
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
   linearVelocity_.pop_back();
   invMass_.pop_back();
   force_.pop_back();
   oldForce_.pop_back();
   shapeID_.pop_back();
   rotation_.pop_back();
   angularVelocity_.pop_back();
   torque_.pop_back();
   oldTorque_.pop_back();
   radiusAtTemperature_.pop_back();
   currentBlock_.pop_back();
   type_.pop_back();
   nextParticle_.pop_back();
   oldContactHistory_.pop_back();
   newContactHistory_.pop_back();
   temperature_.pop_back();
   heatFlux_.pop_back();
   dv_.pop_back();
   dw_.pop_back();
   hydrodynamicForce_.pop_back();
   hydrodynamicTorque_.pop_back();
   oldHydrodynamicForce_.pop_back();
   oldHydrodynamicTorque_.pop_back();
   neighborState_.pop_back();
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
   uid_.reserve(size);
   position_.reserve(size);
   interactionRadius_.reserve(size);
   flags_.reserve(size);
   owner_.reserve(size);
   ghostOwners_.reserve(size);
   linearVelocity_.reserve(size);
   invMass_.reserve(size);
   force_.reserve(size);
   oldForce_.reserve(size);
   shapeID_.reserve(size);
   rotation_.reserve(size);
   angularVelocity_.reserve(size);
   torque_.reserve(size);
   oldTorque_.reserve(size);
   radiusAtTemperature_.reserve(size);
   currentBlock_.reserve(size);
   type_.reserve(size);
   nextParticle_.reserve(size);
   oldContactHistory_.reserve(size);
   newContactHistory_.reserve(size);
   temperature_.reserve(size);
   heatFlux_.reserve(size);
   dv_.reserve(size);
   dw_.reserve(size);
   hydrodynamicForce_.reserve(size);
   hydrodynamicTorque_.reserve(size);
   oldHydrodynamicForce_.reserve(size);
   oldHydrodynamicTorque_.reserve(size);
   neighborState_.reserve(size);
}

inline void ParticleStorage::clear()
{
   uid_.clear();
   position_.clear();
   interactionRadius_.clear();
   flags_.clear();
   owner_.clear();
   ghostOwners_.clear();
   linearVelocity_.clear();
   invMass_.clear();
   force_.clear();
   oldForce_.clear();
   shapeID_.clear();
   rotation_.clear();
   angularVelocity_.clear();
   torque_.clear();
   oldTorque_.clear();
   radiusAtTemperature_.clear();
   currentBlock_.clear();
   type_.clear();
   nextParticle_.clear();
   oldContactHistory_.clear();
   newContactHistory_.clear();
   temperature_.clear();
   heatFlux_.clear();
   dv_.clear();
   dw_.clear();
   hydrodynamicForce_.clear();
   hydrodynamicTorque_.clear();
   oldHydrodynamicForce_.clear();
   oldHydrodynamicTorque_.clear();
   neighborState_.clear();
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
   //WALBERLA_ASSERT_EQUAL( uid_.size(), linearVelocity.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), invMass.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), force.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldForce.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), shapeID.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), rotation.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), angularVelocity.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), torque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldTorque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), radiusAtTemperature.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), currentBlock.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), type.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), nextParticle.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldContactHistory.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), newContactHistory.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), temperature.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), heatFlux.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), dv.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), dw.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), hydrodynamicForce.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), hydrodynamicTorque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldHydrodynamicForce.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), oldHydrodynamicTorque.size() );
   //WALBERLA_ASSERT_EQUAL( uid_.size(), neighborState.size() );
   return uid_.size();
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
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void ParticleStorage::forEachParticle(const bool openmp,
                                             const Selector& selector,
                                             Accessor& acForPS,
                                             Func&& func, Args&&... args) 
{
   static_assert (std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor" );

   WALBERLA_UNUSED(openmp);
   const uint64_t len = size();
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
   walberla::id_t const & operator()(const data::Particle& p) const {return p.getUid();}
};
///Predicate that selects a certain property from a Particle
class SelectParticlePosition
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getPositionRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getPositionRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getPosition();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleInteractionRadius
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getInteractionRadiusRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getInteractionRadiusRef();}
   walberla::real_t const & operator()(const data::Particle& p) const {return p.getInteractionRadius();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleFlags
{
public:
   using return_type = walberla::mesa_pd::data::particle_flags::FlagT;
   walberla::mesa_pd::data::particle_flags::FlagT& operator()(data::Particle& p) const {return p.getFlagsRef();}
   walberla::mesa_pd::data::particle_flags::FlagT& operator()(data::Particle&& p) const {return p.getFlagsRef();}
   walberla::mesa_pd::data::particle_flags::FlagT const & operator()(const data::Particle& p) const {return p.getFlags();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOwner
{
public:
   using return_type = int;
   int& operator()(data::Particle& p) const {return p.getOwnerRef();}
   int& operator()(data::Particle&& p) const {return p.getOwnerRef();}
   int const & operator()(const data::Particle& p) const {return p.getOwner();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleGhostOwners
{
public:
   using return_type = std::unordered_set<walberla::mpi::MPIRank>;
   std::unordered_set<walberla::mpi::MPIRank>& operator()(data::Particle& p) const {return p.getGhostOwnersRef();}
   std::unordered_set<walberla::mpi::MPIRank>& operator()(data::Particle&& p) const {return p.getGhostOwnersRef();}
   std::unordered_set<walberla::mpi::MPIRank> const & operator()(const data::Particle& p) const {return p.getGhostOwners();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleLinearVelocity
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getLinearVelocityRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getLinearVelocityRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getLinearVelocity();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleInvMass
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getInvMassRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getInvMassRef();}
   walberla::real_t const & operator()(const data::Particle& p) const {return p.getInvMass();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getForceRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldForceRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getOldForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleShapeID
{
public:
   using return_type = size_t;
   size_t& operator()(data::Particle& p) const {return p.getShapeIDRef();}
   size_t& operator()(data::Particle&& p) const {return p.getShapeIDRef();}
   size_t const & operator()(const data::Particle& p) const {return p.getShapeID();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleRotation
{
public:
   using return_type = walberla::mesa_pd::Rot3;
   walberla::mesa_pd::Rot3& operator()(data::Particle& p) const {return p.getRotationRef();}
   walberla::mesa_pd::Rot3& operator()(data::Particle&& p) const {return p.getRotationRef();}
   walberla::mesa_pd::Rot3 const & operator()(const data::Particle& p) const {return p.getRotation();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleAngularVelocity
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getAngularVelocityRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getAngularVelocityRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getAngularVelocity();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getTorqueRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldTorqueRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getOldTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleRadiusAtTemperature
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getRadiusAtTemperatureRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getRadiusAtTemperatureRef();}
   walberla::real_t const & operator()(const data::Particle& p) const {return p.getRadiusAtTemperature();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleCurrentBlock
{
public:
   using return_type = blockforest::BlockID;
   blockforest::BlockID& operator()(data::Particle& p) const {return p.getCurrentBlockRef();}
   blockforest::BlockID& operator()(data::Particle&& p) const {return p.getCurrentBlockRef();}
   blockforest::BlockID const & operator()(const data::Particle& p) const {return p.getCurrentBlock();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleType
{
public:
   using return_type = uint_t;
   uint_t& operator()(data::Particle& p) const {return p.getTypeRef();}
   uint_t& operator()(data::Particle&& p) const {return p.getTypeRef();}
   uint_t const & operator()(const data::Particle& p) const {return p.getType();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleNextParticle
{
public:
   using return_type = int;
   int& operator()(data::Particle& p) const {return p.getNextParticleRef();}
   int& operator()(data::Particle&& p) const {return p.getNextParticleRef();}
   int const & operator()(const data::Particle& p) const {return p.getNextParticle();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldContactHistory
{
public:
   using return_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle& p) const {return p.getOldContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle&& p) const {return p.getOldContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & operator()(const data::Particle& p) const {return p.getOldContactHistory();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleNewContactHistory
{
public:
   using return_type = std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle& p) const {return p.getNewContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& operator()(data::Particle&& p) const {return p.getNewContactHistoryRef();}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & operator()(const data::Particle& p) const {return p.getNewContactHistory();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleTemperature
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getTemperatureRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getTemperatureRef();}
   walberla::real_t const & operator()(const data::Particle& p) const {return p.getTemperature();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleHeatFlux
{
public:
   using return_type = walberla::real_t;
   walberla::real_t& operator()(data::Particle& p) const {return p.getHeatFluxRef();}
   walberla::real_t& operator()(data::Particle&& p) const {return p.getHeatFluxRef();}
   walberla::real_t const & operator()(const data::Particle& p) const {return p.getHeatFlux();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleDv
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getDvRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getDvRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getDv();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleDw
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getDwRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getDwRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getDw();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleHydrodynamicForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getHydrodynamicForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getHydrodynamicForceRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getHydrodynamicForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleHydrodynamicTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getHydrodynamicTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getHydrodynamicTorqueRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getHydrodynamicTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldHydrodynamicForce
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldHydrodynamicForceRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldHydrodynamicForceRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getOldHydrodynamicForce();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleOldHydrodynamicTorque
{
public:
   using return_type = walberla::mesa_pd::Vec3;
   walberla::mesa_pd::Vec3& operator()(data::Particle& p) const {return p.getOldHydrodynamicTorqueRef();}
   walberla::mesa_pd::Vec3& operator()(data::Particle&& p) const {return p.getOldHydrodynamicTorqueRef();}
   walberla::mesa_pd::Vec3 const & operator()(const data::Particle& p) const {return p.getOldHydrodynamicTorque();}
};
///Predicate that selects a certain property from a Particle
class SelectParticleNeighborState
{
public:
   using return_type = std::unordered_set<walberla::mpi::MPIRank>;
   std::unordered_set<walberla::mpi::MPIRank>& operator()(data::Particle& p) const {return p.getNeighborStateRef();}
   std::unordered_set<walberla::mpi::MPIRank>& operator()(data::Particle&& p) const {return p.getNeighborStateRef();}
   std::unordered_set<walberla::mpi::MPIRank> const & operator()(const data::Particle& p) const {return p.getNeighborState();}
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla