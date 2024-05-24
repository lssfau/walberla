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
   ~ParticleAccessor() override = default;
   walberla::id_t const & getUid(const size_t p_idx) const {return ps_->getUid(p_idx);}
   walberla::id_t& getUidRef(const size_t p_idx) {return ps_->getUidRef(p_idx);}
   void setUid(const size_t p_idx, walberla::id_t const & v) { ps_->setUid(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getPosition(const size_t p_idx) const {return ps_->getPosition(p_idx);}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t p_idx) {return ps_->getPositionRef(p_idx);}
   void setPosition(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setPosition(p_idx, v);}
   
   walberla::real_t const & getInteractionRadius(const size_t p_idx) const {return ps_->getInteractionRadius(p_idx);}
   walberla::real_t& getInteractionRadiusRef(const size_t p_idx) {return ps_->getInteractionRadiusRef(p_idx);}
   void setInteractionRadius(const size_t p_idx, walberla::real_t const & v) { ps_->setInteractionRadius(p_idx, v);}
   
   walberla::mesa_pd::data::particle_flags::FlagT const & getFlags(const size_t p_idx) const {return ps_->getFlags(p_idx);}
   walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef(const size_t p_idx) {return ps_->getFlagsRef(p_idx);}
   void setFlags(const size_t p_idx, walberla::mesa_pd::data::particle_flags::FlagT const & v) { ps_->setFlags(p_idx, v);}
   
   int const & getOwner(const size_t p_idx) const {return ps_->getOwner(p_idx);}
   int& getOwnerRef(const size_t p_idx) {return ps_->getOwnerRef(p_idx);}
   void setOwner(const size_t p_idx, int const & v) { ps_->setOwner(p_idx, v);}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t p_idx) const {return ps_->getGhostOwners(p_idx);}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t p_idx) {return ps_->getGhostOwnersRef(p_idx);}
   void setGhostOwners(const size_t p_idx, std::unordered_set<walberla::mpi::MPIRank> const & v) { ps_->setGhostOwners(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getLinearVelocity(const size_t p_idx) const {return ps_->getLinearVelocity(p_idx);}
   walberla::mesa_pd::Vec3& getLinearVelocityRef(const size_t p_idx) {return ps_->getLinearVelocityRef(p_idx);}
   void setLinearVelocity(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setLinearVelocity(p_idx, v);}
   
   walberla::real_t const & getInvMass(const size_t p_idx) const {return ps_->getInvMass(p_idx);}
   walberla::real_t& getInvMassRef(const size_t p_idx) {return ps_->getInvMassRef(p_idx);}
   void setInvMass(const size_t p_idx, walberla::real_t const & v) { ps_->setInvMass(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getForce(const size_t p_idx) const {return ps_->getForce(p_idx);}
   walberla::mesa_pd::Vec3& getForceRef(const size_t p_idx) {return ps_->getForceRef(p_idx);}
   void setForce(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setForce(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldForce(const size_t p_idx) const {return ps_->getOldForce(p_idx);}
   walberla::mesa_pd::Vec3& getOldForceRef(const size_t p_idx) {return ps_->getOldForceRef(p_idx);}
   void setOldForce(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldForce(p_idx, v);}
   
   size_t const & getShapeID(const size_t p_idx) const {return ps_->getShapeID(p_idx);}
   size_t& getShapeIDRef(const size_t p_idx) {return ps_->getShapeIDRef(p_idx);}
   void setShapeID(const size_t p_idx, size_t const & v) { ps_->setShapeID(p_idx, v);}
   
   std::shared_ptr<walberla::mesa_pd::data::BaseShape> const & getBaseShape(const size_t p_idx) const {return ps_->getBaseShape(p_idx);}
   std::shared_ptr<walberla::mesa_pd::data::BaseShape>& getBaseShapeRef(const size_t p_idx) {return ps_->getBaseShapeRef(p_idx);}
   void setBaseShape(const size_t p_idx, std::shared_ptr<walberla::mesa_pd::data::BaseShape> const & v) { ps_->setBaseShape(p_idx, v);}
   
   walberla::mesa_pd::Rot3 const & getRotation(const size_t p_idx) const {return ps_->getRotation(p_idx);}
   walberla::mesa_pd::Rot3& getRotationRef(const size_t p_idx) {return ps_->getRotationRef(p_idx);}
   void setRotation(const size_t p_idx, walberla::mesa_pd::Rot3 const & v) { ps_->setRotation(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getAngularVelocity(const size_t p_idx) const {return ps_->getAngularVelocity(p_idx);}
   walberla::mesa_pd::Vec3& getAngularVelocityRef(const size_t p_idx) {return ps_->getAngularVelocityRef(p_idx);}
   void setAngularVelocity(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setAngularVelocity(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getTorque(const size_t p_idx) const {return ps_->getTorque(p_idx);}
   walberla::mesa_pd::Vec3& getTorqueRef(const size_t p_idx) {return ps_->getTorqueRef(p_idx);}
   void setTorque(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setTorque(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldTorque(const size_t p_idx) const {return ps_->getOldTorque(p_idx);}
   walberla::mesa_pd::Vec3& getOldTorqueRef(const size_t p_idx) {return ps_->getOldTorqueRef(p_idx);}
   void setOldTorque(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldTorque(p_idx, v);}
   
   walberla::real_t const & getRadiusAtTemperature(const size_t p_idx) const {return ps_->getRadiusAtTemperature(p_idx);}
   walberla::real_t& getRadiusAtTemperatureRef(const size_t p_idx) {return ps_->getRadiusAtTemperatureRef(p_idx);}
   void setRadiusAtTemperature(const size_t p_idx, walberla::real_t const & v) { ps_->setRadiusAtTemperature(p_idx, v);}
   
   blockforest::BlockID const & getCurrentBlock(const size_t p_idx) const {return ps_->getCurrentBlock(p_idx);}
   blockforest::BlockID& getCurrentBlockRef(const size_t p_idx) {return ps_->getCurrentBlockRef(p_idx);}
   void setCurrentBlock(const size_t p_idx, blockforest::BlockID const & v) { ps_->setCurrentBlock(p_idx, v);}
   
   uint_t const & getType(const size_t p_idx) const {return ps_->getType(p_idx);}
   uint_t& getTypeRef(const size_t p_idx) {return ps_->getTypeRef(p_idx);}
   void setType(const size_t p_idx, uint_t const & v) { ps_->setType(p_idx, v);}
   
   int const & getNextParticle(const size_t p_idx) const {return ps_->getNextParticle(p_idx);}
   int& getNextParticleRef(const size_t p_idx) {return ps_->getNextParticleRef(p_idx);}
   void setNextParticle(const size_t p_idx, int const & v) { ps_->setNextParticle(p_idx, v);}
   
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & getOldContactHistory(const size_t p_idx) const {return ps_->getOldContactHistory(p_idx);}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef(const size_t p_idx) {return ps_->getOldContactHistoryRef(p_idx);}
   void setOldContactHistory(const size_t p_idx, std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & v) { ps_->setOldContactHistory(p_idx, v);}
   
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & getNewContactHistory(const size_t p_idx) const {return ps_->getNewContactHistory(p_idx);}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef(const size_t p_idx) {return ps_->getNewContactHistoryRef(p_idx);}
   void setNewContactHistory(const size_t p_idx, std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & v) { ps_->setNewContactHistory(p_idx, v);}
   
   walberla::real_t const & getTemperature(const size_t p_idx) const {return ps_->getTemperature(p_idx);}
   walberla::real_t& getTemperatureRef(const size_t p_idx) {return ps_->getTemperatureRef(p_idx);}
   void setTemperature(const size_t p_idx, walberla::real_t const & v) { ps_->setTemperature(p_idx, v);}
   
   walberla::real_t const & getHeatFlux(const size_t p_idx) const {return ps_->getHeatFlux(p_idx);}
   walberla::real_t& getHeatFluxRef(const size_t p_idx) {return ps_->getHeatFluxRef(p_idx);}
   void setHeatFlux(const size_t p_idx, walberla::real_t const & v) { ps_->setHeatFlux(p_idx, v);}
   
   uint_t const & getNumContacts(const size_t p_idx) const {return ps_->getNumContacts(p_idx);}
   uint_t& getNumContactsRef(const size_t p_idx) {return ps_->getNumContactsRef(p_idx);}
   void setNumContacts(const size_t p_idx, uint_t const & v) { ps_->setNumContacts(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getDv(const size_t p_idx) const {return ps_->getDv(p_idx);}
   walberla::mesa_pd::Vec3& getDvRef(const size_t p_idx) {return ps_->getDvRef(p_idx);}
   void setDv(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setDv(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getDw(const size_t p_idx) const {return ps_->getDw(p_idx);}
   walberla::mesa_pd::Vec3& getDwRef(const size_t p_idx) {return ps_->getDwRef(p_idx);}
   void setDw(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setDw(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getHydrodynamicForce(const size_t p_idx) const {return ps_->getHydrodynamicForce(p_idx);}
   walberla::mesa_pd::Vec3& getHydrodynamicForceRef(const size_t p_idx) {return ps_->getHydrodynamicForceRef(p_idx);}
   void setHydrodynamicForce(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setHydrodynamicForce(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getHydrodynamicTorque(const size_t p_idx) const {return ps_->getHydrodynamicTorque(p_idx);}
   walberla::mesa_pd::Vec3& getHydrodynamicTorqueRef(const size_t p_idx) {return ps_->getHydrodynamicTorqueRef(p_idx);}
   void setHydrodynamicTorque(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setHydrodynamicTorque(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldHydrodynamicForce(const size_t p_idx) const {return ps_->getOldHydrodynamicForce(p_idx);}
   walberla::mesa_pd::Vec3& getOldHydrodynamicForceRef(const size_t p_idx) {return ps_->getOldHydrodynamicForceRef(p_idx);}
   void setOldHydrodynamicForce(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldHydrodynamicForce(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldHydrodynamicTorque(const size_t p_idx) const {return ps_->getOldHydrodynamicTorque(p_idx);}
   walberla::mesa_pd::Vec3& getOldHydrodynamicTorqueRef(const size_t p_idx) {return ps_->getOldHydrodynamicTorqueRef(p_idx);}
   void setOldHydrodynamicTorque(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldHydrodynamicTorque(p_idx, v);}
   
   walberla::real_t const & getTotalDisplacement(const size_t p_idx) const {return ps_->getTotalDisplacement(p_idx);}
   walberla::real_t& getTotalDisplacementRef(const size_t p_idx) {return ps_->getTotalDisplacementRef(p_idx);}
   void setTotalDisplacement(const size_t p_idx, walberla::real_t const & v) { ps_->setTotalDisplacement(p_idx, v);}
   
   walberla::real_t const & getCollisionForceNorm(const size_t p_idx) const {return ps_->getCollisionForceNorm(p_idx);}
   walberla::real_t& getCollisionForceNormRef(const size_t p_idx) {return ps_->getCollisionForceNormRef(p_idx);}
   void setCollisionForceNorm(const size_t p_idx, walberla::real_t const & v) { ps_->setCollisionForceNorm(p_idx, v);}
   
   walberla::real_t const & getVirtualMass(const size_t p_idx) const {return ps_->getVirtualMass(p_idx);}
   walberla::real_t& getVirtualMassRef(const size_t p_idx) {return ps_->getVirtualMassRef(p_idx);}
   void setVirtualMass(const size_t p_idx, walberla::real_t const & v) { ps_->setVirtualMass(p_idx, v);}
   
   walberla::real_t const & getInvMassIncludingVirtual(const size_t p_idx) const {return ps_->getInvMassIncludingVirtual(p_idx);}
   walberla::real_t& getInvMassIncludingVirtualRef(const size_t p_idx) {return ps_->getInvMassIncludingVirtualRef(p_idx);}
   void setInvMassIncludingVirtual(const size_t p_idx, walberla::real_t const & v) { ps_->setInvMassIncludingVirtual(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldLinearAcceleration(const size_t p_idx) const {return ps_->getOldLinearAcceleration(p_idx);}
   walberla::mesa_pd::Vec3& getOldLinearAccelerationRef(const size_t p_idx) {return ps_->getOldLinearAccelerationRef(p_idx);}
   void setOldLinearAcceleration(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldLinearAcceleration(p_idx, v);}
   
   walberla::mesa_pd::Mat3 const & getInvInertiaBF(const size_t p_idx) const {return ps_->getInvInertiaBF(p_idx);}
   walberla::mesa_pd::Mat3& getInvInertiaBFRef(const size_t p_idx) {return ps_->getInvInertiaBFRef(p_idx);}
   void setInvInertiaBF(const size_t p_idx, walberla::mesa_pd::Mat3 const & v) { ps_->setInvInertiaBF(p_idx, v);}
   
   walberla::mesa_pd::Mat3 const & getVirtualInertiaBF(const size_t p_idx) const {return ps_->getVirtualInertiaBF(p_idx);}
   walberla::mesa_pd::Mat3& getVirtualInertiaBFRef(const size_t p_idx) {return ps_->getVirtualInertiaBFRef(p_idx);}
   void setVirtualInertiaBF(const size_t p_idx, walberla::mesa_pd::Mat3 const & v) { ps_->setVirtualInertiaBF(p_idx, v);}
   
   walberla::mesa_pd::Mat3 const & getInvInertiaBFIncludingVirtual(const size_t p_idx) const {return ps_->getInvInertiaBFIncludingVirtual(p_idx);}
   walberla::mesa_pd::Mat3& getInvInertiaBFIncludingVirtualRef(const size_t p_idx) {return ps_->getInvInertiaBFIncludingVirtualRef(p_idx);}
   void setInvInertiaBFIncludingVirtual(const size_t p_idx, walberla::mesa_pd::Mat3 const & v) { ps_->setInvInertiaBFIncludingVirtual(p_idx, v);}
   
   walberla::mesa_pd::Vec3 const & getOldAngularAcceleration(const size_t p_idx) const {return ps_->getOldAngularAcceleration(p_idx);}
   walberla::mesa_pd::Vec3& getOldAngularAccelerationRef(const size_t p_idx) {return ps_->getOldAngularAccelerationRef(p_idx);}
   void setOldAngularAcceleration(const size_t p_idx, walberla::mesa_pd::Vec3 const & v) { ps_->setOldAngularAcceleration(p_idx, v);}
   
   int64_t const & getClusterID(const size_t p_idx) const {return ps_->getClusterID(p_idx);}
   int64_t& getClusterIDRef(const size_t p_idx) {return ps_->getClusterIDRef(p_idx);}
   void setClusterID(const size_t p_idx, int64_t const & v) { ps_->setClusterID(p_idx, v);}
   
   int64_t const & getSegmentID(const size_t p_idx) const {return ps_->getSegmentID(p_idx);}
   int64_t& getSegmentIDRef(const size_t p_idx) {return ps_->getSegmentIDRef(p_idx);}
   void setSegmentID(const size_t p_idx, int64_t const & v) { ps_->setSegmentID(p_idx, v);}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getNeighborState(const size_t p_idx) const {return ps_->getNeighborState(p_idx);}
   std::unordered_set<walberla::mpi::MPIRank>& getNeighborStateRef(const size_t p_idx) {return ps_->getNeighborStateRef(p_idx);}
   void setNeighborState(const size_t p_idx, std::unordered_set<walberla::mpi::MPIRank> const & v) { ps_->setNeighborState(p_idx, v);}
   

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
   ~SingleParticleAccessor() override = default;
   walberla::id_t const & getUid(const size_t /*p_idx*/) const {return uid_;}
   void setUid(const size_t /*p_idx*/, walberla::id_t const & v) { uid_ = v;}
   walberla::id_t& getUidRef(const size_t /*p_idx*/) {return uid_;}
   
   walberla::mesa_pd::Vec3 const & getPosition(const size_t /*p_idx*/) const {return position_;}
   void setPosition(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { position_ = v;}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t /*p_idx*/) {return position_;}
   
   walberla::real_t const & getInteractionRadius(const size_t /*p_idx*/) const {return interactionRadius_;}
   void setInteractionRadius(const size_t /*p_idx*/, walberla::real_t const & v) { interactionRadius_ = v;}
   walberla::real_t& getInteractionRadiusRef(const size_t /*p_idx*/) {return interactionRadius_;}
   
   walberla::mesa_pd::data::particle_flags::FlagT const & getFlags(const size_t /*p_idx*/) const {return flags_;}
   void setFlags(const size_t /*p_idx*/, walberla::mesa_pd::data::particle_flags::FlagT const & v) { flags_ = v;}
   walberla::mesa_pd::data::particle_flags::FlagT& getFlagsRef(const size_t /*p_idx*/) {return flags_;}
   
   int const & getOwner(const size_t /*p_idx*/) const {return owner_;}
   void setOwner(const size_t /*p_idx*/, int const & v) { owner_ = v;}
   int& getOwnerRef(const size_t /*p_idx*/) {return owner_;}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t /*p_idx*/) const {return ghostOwners_;}
   void setGhostOwners(const size_t /*p_idx*/, std::unordered_set<walberla::mpi::MPIRank> const & v) { ghostOwners_ = v;}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t /*p_idx*/) {return ghostOwners_;}
   
   walberla::mesa_pd::Vec3 const & getLinearVelocity(const size_t /*p_idx*/) const {return linearVelocity_;}
   void setLinearVelocity(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { linearVelocity_ = v;}
   walberla::mesa_pd::Vec3& getLinearVelocityRef(const size_t /*p_idx*/) {return linearVelocity_;}
   
   walberla::real_t const & getInvMass(const size_t /*p_idx*/) const {return invMass_;}
   void setInvMass(const size_t /*p_idx*/, walberla::real_t const & v) { invMass_ = v;}
   walberla::real_t& getInvMassRef(const size_t /*p_idx*/) {return invMass_;}
   
   walberla::mesa_pd::Vec3 const & getForce(const size_t /*p_idx*/) const {return force_;}
   void setForce(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { force_ = v;}
   walberla::mesa_pd::Vec3& getForceRef(const size_t /*p_idx*/) {return force_;}
   
   walberla::mesa_pd::Vec3 const & getOldForce(const size_t /*p_idx*/) const {return oldForce_;}
   void setOldForce(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldForce_ = v;}
   walberla::mesa_pd::Vec3& getOldForceRef(const size_t /*p_idx*/) {return oldForce_;}
   
   size_t const & getShapeID(const size_t /*p_idx*/) const {return shapeID_;}
   void setShapeID(const size_t /*p_idx*/, size_t const & v) { shapeID_ = v;}
   size_t& getShapeIDRef(const size_t /*p_idx*/) {return shapeID_;}
   
   std::shared_ptr<walberla::mesa_pd::data::BaseShape> const & getBaseShape(const size_t /*p_idx*/) const {return baseShape_;}
   void setBaseShape(const size_t /*p_idx*/, std::shared_ptr<walberla::mesa_pd::data::BaseShape> const & v) { baseShape_ = v;}
   std::shared_ptr<walberla::mesa_pd::data::BaseShape>& getBaseShapeRef(const size_t /*p_idx*/) {return baseShape_;}
   
   walberla::mesa_pd::Rot3 const & getRotation(const size_t /*p_idx*/) const {return rotation_;}
   void setRotation(const size_t /*p_idx*/, walberla::mesa_pd::Rot3 const & v) { rotation_ = v;}
   walberla::mesa_pd::Rot3& getRotationRef(const size_t /*p_idx*/) {return rotation_;}
   
   walberla::mesa_pd::Vec3 const & getAngularVelocity(const size_t /*p_idx*/) const {return angularVelocity_;}
   void setAngularVelocity(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { angularVelocity_ = v;}
   walberla::mesa_pd::Vec3& getAngularVelocityRef(const size_t /*p_idx*/) {return angularVelocity_;}
   
   walberla::mesa_pd::Vec3 const & getTorque(const size_t /*p_idx*/) const {return torque_;}
   void setTorque(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { torque_ = v;}
   walberla::mesa_pd::Vec3& getTorqueRef(const size_t /*p_idx*/) {return torque_;}
   
   walberla::mesa_pd::Vec3 const & getOldTorque(const size_t /*p_idx*/) const {return oldTorque_;}
   void setOldTorque(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldTorque_ = v;}
   walberla::mesa_pd::Vec3& getOldTorqueRef(const size_t /*p_idx*/) {return oldTorque_;}
   
   walberla::real_t const & getRadiusAtTemperature(const size_t /*p_idx*/) const {return radiusAtTemperature_;}
   void setRadiusAtTemperature(const size_t /*p_idx*/, walberla::real_t const & v) { radiusAtTemperature_ = v;}
   walberla::real_t& getRadiusAtTemperatureRef(const size_t /*p_idx*/) {return radiusAtTemperature_;}
   
   blockforest::BlockID const & getCurrentBlock(const size_t /*p_idx*/) const {return currentBlock_;}
   void setCurrentBlock(const size_t /*p_idx*/, blockforest::BlockID const & v) { currentBlock_ = v;}
   blockforest::BlockID& getCurrentBlockRef(const size_t /*p_idx*/) {return currentBlock_;}
   
   uint_t const & getType(const size_t /*p_idx*/) const {return type_;}
   void setType(const size_t /*p_idx*/, uint_t const & v) { type_ = v;}
   uint_t& getTypeRef(const size_t /*p_idx*/) {return type_;}
   
   int const & getNextParticle(const size_t /*p_idx*/) const {return nextParticle_;}
   void setNextParticle(const size_t /*p_idx*/, int const & v) { nextParticle_ = v;}
   int& getNextParticleRef(const size_t /*p_idx*/) {return nextParticle_;}
   
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & getOldContactHistory(const size_t /*p_idx*/) const {return oldContactHistory_;}
   void setOldContactHistory(const size_t /*p_idx*/, std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & v) { oldContactHistory_ = v;}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getOldContactHistoryRef(const size_t /*p_idx*/) {return oldContactHistory_;}
   
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & getNewContactHistory(const size_t /*p_idx*/) const {return newContactHistory_;}
   void setNewContactHistory(const size_t /*p_idx*/, std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> const & v) { newContactHistory_ = v;}
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>& getNewContactHistoryRef(const size_t /*p_idx*/) {return newContactHistory_;}
   
   walberla::real_t const & getTemperature(const size_t /*p_idx*/) const {return temperature_;}
   void setTemperature(const size_t /*p_idx*/, walberla::real_t const & v) { temperature_ = v;}
   walberla::real_t& getTemperatureRef(const size_t /*p_idx*/) {return temperature_;}
   
   walberla::real_t const & getHeatFlux(const size_t /*p_idx*/) const {return heatFlux_;}
   void setHeatFlux(const size_t /*p_idx*/, walberla::real_t const & v) { heatFlux_ = v;}
   walberla::real_t& getHeatFluxRef(const size_t /*p_idx*/) {return heatFlux_;}
   
   uint_t const & getNumContacts(const size_t /*p_idx*/) const {return numContacts_;}
   void setNumContacts(const size_t /*p_idx*/, uint_t const & v) { numContacts_ = v;}
   uint_t& getNumContactsRef(const size_t /*p_idx*/) {return numContacts_;}
   
   walberla::mesa_pd::Vec3 const & getDv(const size_t /*p_idx*/) const {return dv_;}
   void setDv(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { dv_ = v;}
   walberla::mesa_pd::Vec3& getDvRef(const size_t /*p_idx*/) {return dv_;}
   
   walberla::mesa_pd::Vec3 const & getDw(const size_t /*p_idx*/) const {return dw_;}
   void setDw(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { dw_ = v;}
   walberla::mesa_pd::Vec3& getDwRef(const size_t /*p_idx*/) {return dw_;}
   
   walberla::mesa_pd::Vec3 const & getHydrodynamicForce(const size_t /*p_idx*/) const {return hydrodynamicForce_;}
   void setHydrodynamicForce(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { hydrodynamicForce_ = v;}
   walberla::mesa_pd::Vec3& getHydrodynamicForceRef(const size_t /*p_idx*/) {return hydrodynamicForce_;}
   
   walberla::mesa_pd::Vec3 const & getHydrodynamicTorque(const size_t /*p_idx*/) const {return hydrodynamicTorque_;}
   void setHydrodynamicTorque(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { hydrodynamicTorque_ = v;}
   walberla::mesa_pd::Vec3& getHydrodynamicTorqueRef(const size_t /*p_idx*/) {return hydrodynamicTorque_;}
   
   walberla::mesa_pd::Vec3 const & getOldHydrodynamicForce(const size_t /*p_idx*/) const {return oldHydrodynamicForce_;}
   void setOldHydrodynamicForce(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldHydrodynamicForce_ = v;}
   walberla::mesa_pd::Vec3& getOldHydrodynamicForceRef(const size_t /*p_idx*/) {return oldHydrodynamicForce_;}
   
   walberla::mesa_pd::Vec3 const & getOldHydrodynamicTorque(const size_t /*p_idx*/) const {return oldHydrodynamicTorque_;}
   void setOldHydrodynamicTorque(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldHydrodynamicTorque_ = v;}
   walberla::mesa_pd::Vec3& getOldHydrodynamicTorqueRef(const size_t /*p_idx*/) {return oldHydrodynamicTorque_;}
   
   walberla::real_t const & getTotalDisplacement(const size_t /*p_idx*/) const {return totalDisplacement_;}
   void setTotalDisplacement(const size_t /*p_idx*/, walberla::real_t const & v) { totalDisplacement_ = v;}
   walberla::real_t& getTotalDisplacementRef(const size_t /*p_idx*/) {return totalDisplacement_;}
   
   walberla::real_t const & getCollisionForceNorm(const size_t /*p_idx*/) const {return collisionForceNorm_;}
   void setCollisionForceNorm(const size_t /*p_idx*/, walberla::real_t const & v) { collisionForceNorm_ = v;}
   walberla::real_t& getCollisionForceNormRef(const size_t /*p_idx*/) {return collisionForceNorm_;}
   
   walberla::real_t const & getVirtualMass(const size_t /*p_idx*/) const {return virtualMass_;}
   void setVirtualMass(const size_t /*p_idx*/, walberla::real_t const & v) { virtualMass_ = v;}
   walberla::real_t& getVirtualMassRef(const size_t /*p_idx*/) {return virtualMass_;}
   
   walberla::real_t const & getInvMassIncludingVirtual(const size_t /*p_idx*/) const {return invMassIncludingVirtual_;}
   void setInvMassIncludingVirtual(const size_t /*p_idx*/, walberla::real_t const & v) { invMassIncludingVirtual_ = v;}
   walberla::real_t& getInvMassIncludingVirtualRef(const size_t /*p_idx*/) {return invMassIncludingVirtual_;}
   
   walberla::mesa_pd::Vec3 const & getOldLinearAcceleration(const size_t /*p_idx*/) const {return oldLinearAcceleration_;}
   void setOldLinearAcceleration(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldLinearAcceleration_ = v;}
   walberla::mesa_pd::Vec3& getOldLinearAccelerationRef(const size_t /*p_idx*/) {return oldLinearAcceleration_;}
   
   walberla::mesa_pd::Mat3 const & getInvInertiaBF(const size_t /*p_idx*/) const {return invInertiaBF_;}
   void setInvInertiaBF(const size_t /*p_idx*/, walberla::mesa_pd::Mat3 const & v) { invInertiaBF_ = v;}
   walberla::mesa_pd::Mat3& getInvInertiaBFRef(const size_t /*p_idx*/) {return invInertiaBF_;}
   
   walberla::mesa_pd::Mat3 const & getVirtualInertiaBF(const size_t /*p_idx*/) const {return virtualInertiaBF_;}
   void setVirtualInertiaBF(const size_t /*p_idx*/, walberla::mesa_pd::Mat3 const & v) { virtualInertiaBF_ = v;}
   walberla::mesa_pd::Mat3& getVirtualInertiaBFRef(const size_t /*p_idx*/) {return virtualInertiaBF_;}
   
   walberla::mesa_pd::Mat3 const & getInvInertiaBFIncludingVirtual(const size_t /*p_idx*/) const {return invInertiaBFIncludingVirtual_;}
   void setInvInertiaBFIncludingVirtual(const size_t /*p_idx*/, walberla::mesa_pd::Mat3 const & v) { invInertiaBFIncludingVirtual_ = v;}
   walberla::mesa_pd::Mat3& getInvInertiaBFIncludingVirtualRef(const size_t /*p_idx*/) {return invInertiaBFIncludingVirtual_;}
   
   walberla::mesa_pd::Vec3 const & getOldAngularAcceleration(const size_t /*p_idx*/) const {return oldAngularAcceleration_;}
   void setOldAngularAcceleration(const size_t /*p_idx*/, walberla::mesa_pd::Vec3 const & v) { oldAngularAcceleration_ = v;}
   walberla::mesa_pd::Vec3& getOldAngularAccelerationRef(const size_t /*p_idx*/) {return oldAngularAcceleration_;}
   
   int64_t const & getClusterID(const size_t /*p_idx*/) const {return clusterID_;}
   void setClusterID(const size_t /*p_idx*/, int64_t const & v) { clusterID_ = v;}
   int64_t& getClusterIDRef(const size_t /*p_idx*/) {return clusterID_;}
   
   int64_t const & getSegmentID(const size_t /*p_idx*/) const {return segmentID_;}
   void setSegmentID(const size_t /*p_idx*/, int64_t const & v) { segmentID_ = v;}
   int64_t& getSegmentIDRef(const size_t /*p_idx*/) {return segmentID_;}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getNeighborState(const size_t /*p_idx*/) const {return neighborState_;}
   void setNeighborState(const size_t /*p_idx*/, std::unordered_set<walberla::mpi::MPIRank> const & v) { neighborState_ = v;}
   std::unordered_set<walberla::mpi::MPIRank>& getNeighborStateRef(const size_t /*p_idx*/) {return neighborState_;}
   

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
   std::unordered_set<walberla::mpi::MPIRank> ghostOwners_;
   walberla::mesa_pd::Vec3 linearVelocity_;
   walberla::real_t invMass_;
   walberla::mesa_pd::Vec3 force_;
   walberla::mesa_pd::Vec3 oldForce_;
   size_t shapeID_;
   std::shared_ptr<walberla::mesa_pd::data::BaseShape> baseShape_;
   walberla::mesa_pd::Rot3 rotation_;
   walberla::mesa_pd::Vec3 angularVelocity_;
   walberla::mesa_pd::Vec3 torque_;
   walberla::mesa_pd::Vec3 oldTorque_;
   walberla::real_t radiusAtTemperature_;
   blockforest::BlockID currentBlock_;
   uint_t type_;
   int nextParticle_;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> oldContactHistory_;
   std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> newContactHistory_;
   walberla::real_t temperature_;
   walberla::real_t heatFlux_;
   uint_t numContacts_;
   walberla::mesa_pd::Vec3 dv_;
   walberla::mesa_pd::Vec3 dw_;
   walberla::mesa_pd::Vec3 hydrodynamicForce_;
   walberla::mesa_pd::Vec3 hydrodynamicTorque_;
   walberla::mesa_pd::Vec3 oldHydrodynamicForce_;
   walberla::mesa_pd::Vec3 oldHydrodynamicTorque_;
   walberla::real_t totalDisplacement_;
   walberla::real_t collisionForceNorm_;
   walberla::real_t virtualMass_;
   walberla::real_t invMassIncludingVirtual_;
   walberla::mesa_pd::Vec3 oldLinearAcceleration_;
   walberla::mesa_pd::Mat3 invInertiaBF_;
   walberla::mesa_pd::Mat3 virtualInertiaBF_;
   walberla::mesa_pd::Mat3 invInertiaBFIncludingVirtual_;
   walberla::mesa_pd::Vec3 oldAngularAcceleration_;
   int64_t clusterID_;
   int64_t segmentID_;
   std::unordered_set<walberla::mpi::MPIRank> neighborState_;
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla