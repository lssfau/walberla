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
//! \file ParticleAccessor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <core/UniqueID.h>

#include <limits>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * @brief Basic ParticleAccessor for the ParticleStorage
 *
 * Provides get, set and getRef for all members of the ParticleStorage.
 * Can be used as a basis class for a more advanced ParticleAccessor.
 */
class ParticleAccessor : public IAccessor
{
public:
   ParticleAccessor(const std::shared_ptr<data::ParticleStorage>& ps) : ps_(ps) {}
   virtual ~ParticleAccessor() = default;
   const walberla::id_t& getUid(const size_t p_idx) const {return ps_->getUid(p_idx);}
   walberla::id_t& getUidRef(const size_t p_idx) {return ps_->getUidRef(p_idx);}
   void setUid(const size_t p_idx, const walberla::id_t& v) { ps_->setUid(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const {return ps_->getPosition(p_idx);}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t p_idx) {return ps_->getPositionRef(p_idx);}
   void setPosition(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setPosition(p_idx, v);}
   
   const walberla::real_t& getInteractionRadius(const size_t p_idx) const {return ps_->getInteractionRadius(p_idx);}
   walberla::real_t& getInteractionRadiusRef(const size_t p_idx) {return ps_->getInteractionRadiusRef(p_idx);}
   void setInteractionRadius(const size_t p_idx, const walberla::real_t& v) { ps_->setInteractionRadius(p_idx, v);}
   
   const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t p_idx) const {return ps_->getFlags(p_idx);}
   walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef(const size_t p_idx) {return ps_->getFlagsRef(p_idx);}
   void setFlags(const size_t p_idx, const walberla::mesa_pd::data::particle_flags::FlagT& v) { ps_->setFlags(p_idx, v);}
   
   const int& getOwner(const size_t p_idx) const {return ps_->getOwner(p_idx);}
   int& getOwnerRef(const size_t p_idx) {return ps_->getOwnerRef(p_idx);}
   void setOwner(const size_t p_idx, const int& v) { ps_->setOwner(p_idx, v);}
   
   const std::vector<int>& getGhostOwners(const size_t p_idx) const {return ps_->getGhostOwners(p_idx);}
   std::vector<int>& getGhostOwnersRef(const size_t p_idx) {return ps_->getGhostOwnersRef(p_idx);}
   void setGhostOwners(const size_t p_idx, const std::vector<int>& v) { ps_->setGhostOwners(p_idx, v);}
   
   const size_t& getShapeID(const size_t p_idx) const {return ps_->getShapeID(p_idx);}
   size_t& getShapeIDRef(const size_t p_idx) {return ps_->getShapeIDRef(p_idx);}
   void setShapeID(const size_t p_idx, const size_t& v) { ps_->setShapeID(p_idx, v);}
   
   const walberla::mesa_pd::Rot3& getRotation(const size_t p_idx) const {return ps_->getRotation(p_idx);}
   walberla::mesa_pd::Rot3& getRotationRef(const size_t p_idx) {return ps_->getRotationRef(p_idx);}
   void setRotation(const size_t p_idx, const walberla::mesa_pd::Rot3& v) { ps_->setRotation(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t p_idx) const {return ps_->getAngularVelocity(p_idx);}
   walberla::mesa_pd::Vec3& getAngularVelocityRef(const size_t p_idx) {return ps_->getAngularVelocityRef(p_idx);}
   void setAngularVelocity(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setAngularVelocity(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getTorque(const size_t p_idx) const {return ps_->getTorque(p_idx);}
   walberla::mesa_pd::Vec3& getTorqueRef(const size_t p_idx) {return ps_->getTorqueRef(p_idx);}
   void setTorque(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setTorque(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t p_idx) const {return ps_->getLinearVelocity(p_idx);}
   walberla::mesa_pd::Vec3& getLinearVelocityRef(const size_t p_idx) {return ps_->getLinearVelocityRef(p_idx);}
   void setLinearVelocity(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setLinearVelocity(p_idx, v);}
   
   const walberla::real_t& getInvMass(const size_t p_idx) const {return ps_->getInvMass(p_idx);}
   walberla::real_t& getInvMassRef(const size_t p_idx) {return ps_->getInvMassRef(p_idx);}
   void setInvMass(const size_t p_idx, const walberla::real_t& v) { ps_->setInvMass(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getForce(const size_t p_idx) const {return ps_->getForce(p_idx);}
   walberla::mesa_pd::Vec3& getForceRef(const size_t p_idx) {return ps_->getForceRef(p_idx);}
   void setForce(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setForce(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getOldForce(const size_t p_idx) const {return ps_->getOldForce(p_idx);}
   walberla::mesa_pd::Vec3& getOldForceRef(const size_t p_idx) {return ps_->getOldForceRef(p_idx);}
   void setOldForce(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setOldForce(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getOldTorque(const size_t p_idx) const {return ps_->getOldTorque(p_idx);}
   walberla::mesa_pd::Vec3& getOldTorqueRef(const size_t p_idx) {return ps_->getOldTorqueRef(p_idx);}
   void setOldTorque(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setOldTorque(p_idx, v);}
   
   const uint_t& getType(const size_t p_idx) const {return ps_->getType(p_idx);}
   uint_t& getTypeRef(const size_t p_idx) {return ps_->getTypeRef(p_idx);}
   void setType(const size_t p_idx, const uint_t& v) { ps_->setType(p_idx, v);}
   
   const int& getNextParticle(const size_t p_idx) const {return ps_->getNextParticle(p_idx);}
   int& getNextParticleRef(const size_t p_idx) {return ps_->getNextParticleRef(p_idx);}
   void setNextParticle(const size_t p_idx, const int& v) { ps_->setNextParticle(p_idx, v);}
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistory(const size_t p_idx) const {return ps_->getOldContactHistory(p_idx);}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef(const size_t p_idx) {return ps_->getOldContactHistoryRef(p_idx);}
   void setOldContactHistory(const size_t p_idx, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { ps_->setOldContactHistory(p_idx, v);}
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistory(const size_t p_idx) const {return ps_->getNewContactHistory(p_idx);}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef(const size_t p_idx) {return ps_->getNewContactHistoryRef(p_idx);}
   void setNewContactHistory(const size_t p_idx, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { ps_->setNewContactHistory(p_idx, v);}
   
   const walberla::real_t& getTemperature(const size_t p_idx) const {return ps_->getTemperature(p_idx);}
   walberla::real_t& getTemperatureRef(const size_t p_idx) {return ps_->getTemperatureRef(p_idx);}
   void setTemperature(const size_t p_idx, const walberla::real_t& v) { ps_->setTemperature(p_idx, v);}
   
   const walberla::real_t& getHeatFlux(const size_t p_idx) const {return ps_->getHeatFlux(p_idx);}
   walberla::real_t& getHeatFluxRef(const size_t p_idx) {return ps_->getHeatFluxRef(p_idx);}
   void setHeatFlux(const size_t p_idx, const walberla::real_t& v) { ps_->setHeatFlux(p_idx, v);}
   

   id_t getInvalidUid() const {return UniqueID<data::Particle>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of particle specified by uid.
   * @param uid unique id of the particle to be looked up
   * @return the index of the particle or std::numeric_limits<size_t>::max() if the particle is not found
   */
   size_t uidToIdx(const id_t& uid) const {auto it = ps_->find(uid); return it != ps_->end() ? it.getIdx() : std::numeric_limits<size_t>::max();}
   size_t size() const { return ps_->size(); }

   inline size_t create(const id_t& uid);
   inline size_t erase(const size_t& idx);
   inline size_t find(const id_t& uid);
protected:
   std::shared_ptr<data::ParticleStorage> ps_;
};

inline size_t ParticleAccessor::create(const id_t& uid)
{
   auto it = ps_->create(uid);
   return it.getIdx();
}
inline size_t ParticleAccessor::erase(const size_t& idx)
{
   data::ParticleStorage::iterator it(ps_.get(), idx);
   it = ps_->erase(it);
   return it.getIdx();
}
inline size_t ParticleAccessor::find(const id_t& uid)
{
   auto it = ps_->find(uid);
   return it.getIdx();
}

/**
 * @brief Basic ParticleAccessor which emulates a single particle in a ParticleStorage
 * without actually needing a ParticleStorage. This class is used mainly for testing purposes.
 *
 * Provides get, set and getRef.
 */
class SingleParticleAccessor : public IAccessor
{
public:
   virtual ~SingleParticleAccessor() = default;
   const walberla::id_t& getUid(const size_t /*p_idx*/) const {return uid_;}
   void setUid(const size_t /*p_idx*/, const walberla::id_t& v) { uid_ = v;}
   walberla::id_t& getUidRef(const size_t /*p_idx*/) {return uid_;}
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t /*p_idx*/) const {return position_;}
   void setPosition(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { position_ = v;}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t /*p_idx*/) {return position_;}
   
   const walberla::real_t& getInteractionRadius(const size_t /*p_idx*/) const {return interactionRadius_;}
   void setInteractionRadius(const size_t /*p_idx*/, const walberla::real_t& v) { interactionRadius_ = v;}
   walberla::real_t& getInteractionRadiusRef(const size_t /*p_idx*/) {return interactionRadius_;}
   
   const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t /*p_idx*/) const {return flags_;}
   void setFlags(const size_t /*p_idx*/, const walberla::mesa_pd::data::particle_flags::FlagT& v) { flags_ = v;}
   walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef(const size_t /*p_idx*/) {return flags_;}
   
   const int& getOwner(const size_t /*p_idx*/) const {return owner_;}
   void setOwner(const size_t /*p_idx*/, const int& v) { owner_ = v;}
   int& getOwnerRef(const size_t /*p_idx*/) {return owner_;}
   
   const std::vector<int>& getGhostOwners(const size_t /*p_idx*/) const {return ghostOwners_;}
   void setGhostOwners(const size_t /*p_idx*/, const std::vector<int>& v) { ghostOwners_ = v;}
   std::vector<int>& getGhostOwnersRef(const size_t /*p_idx*/) {return ghostOwners_;}
   
   const size_t& getShapeID(const size_t /*p_idx*/) const {return shapeID_;}
   void setShapeID(const size_t /*p_idx*/, const size_t& v) { shapeID_ = v;}
   size_t& getShapeIDRef(const size_t /*p_idx*/) {return shapeID_;}
   
   const walberla::mesa_pd::Rot3& getRotation(const size_t /*p_idx*/) const {return rotation_;}
   void setRotation(const size_t /*p_idx*/, const walberla::mesa_pd::Rot3& v) { rotation_ = v;}
   walberla::mesa_pd::Rot3& getRotationRef(const size_t /*p_idx*/) {return rotation_;}
   
   const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t /*p_idx*/) const {return angularVelocity_;}
   void setAngularVelocity(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { angularVelocity_ = v;}
   walberla::mesa_pd::Vec3& getAngularVelocityRef(const size_t /*p_idx*/) {return angularVelocity_;}
   
   const walberla::mesa_pd::Vec3& getTorque(const size_t /*p_idx*/) const {return torque_;}
   void setTorque(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { torque_ = v;}
   walberla::mesa_pd::Vec3& getTorqueRef(const size_t /*p_idx*/) {return torque_;}
   
   const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t /*p_idx*/) const {return linearVelocity_;}
   void setLinearVelocity(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { linearVelocity_ = v;}
   walberla::mesa_pd::Vec3& getLinearVelocityRef(const size_t /*p_idx*/) {return linearVelocity_;}
   
   const walberla::real_t& getInvMass(const size_t /*p_idx*/) const {return invMass_;}
   void setInvMass(const size_t /*p_idx*/, const walberla::real_t& v) { invMass_ = v;}
   walberla::real_t& getInvMassRef(const size_t /*p_idx*/) {return invMass_;}
   
   const walberla::mesa_pd::Vec3& getForce(const size_t /*p_idx*/) const {return force_;}
   void setForce(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { force_ = v;}
   walberla::mesa_pd::Vec3& getForceRef(const size_t /*p_idx*/) {return force_;}
   
   const walberla::mesa_pd::Vec3& getOldForce(const size_t /*p_idx*/) const {return oldForce_;}
   void setOldForce(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { oldForce_ = v;}
   walberla::mesa_pd::Vec3& getOldForceRef(const size_t /*p_idx*/) {return oldForce_;}
   
   const walberla::mesa_pd::Vec3& getOldTorque(const size_t /*p_idx*/) const {return oldTorque_;}
   void setOldTorque(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { oldTorque_ = v;}
   walberla::mesa_pd::Vec3& getOldTorqueRef(const size_t /*p_idx*/) {return oldTorque_;}
   
   const uint_t& getType(const size_t /*p_idx*/) const {return type_;}
   void setType(const size_t /*p_idx*/, const uint_t& v) { type_ = v;}
   uint_t& getTypeRef(const size_t /*p_idx*/) {return type_;}
   
   const int& getNextParticle(const size_t /*p_idx*/) const {return nextParticle_;}
   void setNextParticle(const size_t /*p_idx*/, const int& v) { nextParticle_ = v;}
   int& getNextParticleRef(const size_t /*p_idx*/) {return nextParticle_;}
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistory(const size_t /*p_idx*/) const {return oldContactHistory_;}
   void setOldContactHistory(const size_t /*p_idx*/, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { oldContactHistory_ = v;}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef(const size_t /*p_idx*/) {return oldContactHistory_;}
   
   const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistory(const size_t /*p_idx*/) const {return newContactHistory_;}
   void setNewContactHistory(const size_t /*p_idx*/, const std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& v) { newContactHistory_ = v;}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef(const size_t /*p_idx*/) {return newContactHistory_;}
   
   const walberla::real_t& getTemperature(const size_t /*p_idx*/) const {return temperature_;}
   void setTemperature(const size_t /*p_idx*/, const walberla::real_t& v) { temperature_ = v;}
   walberla::real_t& getTemperatureRef(const size_t /*p_idx*/) {return temperature_;}
   
   const walberla::real_t& getHeatFlux(const size_t /*p_idx*/) const {return heatFlux_;}
   void setHeatFlux(const size_t /*p_idx*/, const walberla::real_t& v) { heatFlux_ = v;}
   walberla::real_t& getHeatFluxRef(const size_t /*p_idx*/) {return heatFlux_;}
   

   id_t getInvalidUid() const {return UniqueID<data::Particle>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of particle specified by uid.
   * @param uid unique id of the particle to be looked up
   * @return the index of the particle or std::numeric_limits<size_t>::max() if the particle is not found
   */
   size_t uidToIdx(const id_t& uid) const {return uid == uid_ ? 0 : std::numeric_limits<size_t>::max();}
   size_t size() const { return 1; }
private:
   walberla::id_t uid_;
   walberla::mesa_pd::Vec3 position_;
   walberla::real_t interactionRadius_;
   walberla::mesa_pd::data::particle_flags::FlagT flags_;
   int owner_;
   std::vector<int> ghostOwners_;
   size_t shapeID_;
   walberla::mesa_pd::Rot3 rotation_;
   walberla::mesa_pd::Vec3 angularVelocity_;
   walberla::mesa_pd::Vec3 torque_;
   walberla::mesa_pd::Vec3 linearVelocity_;
   walberla::real_t invMass_;
   walberla::mesa_pd::Vec3 force_;
   walberla::mesa_pd::Vec3 oldForce_;
   walberla::mesa_pd::Vec3 oldTorque_;
   uint_t type_;
   int nextParticle_;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> oldContactHistory_;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> newContactHistory_;
   walberla::real_t temperature_;
   walberla::real_t heatFlux_;
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla