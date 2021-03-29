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
   const walberla::id_t& getUid(const size_t p_idx) const {return ps_->getUid(p_idx);}
   walberla::id_t& getUidRef(const size_t p_idx) {return ps_->getUidRef(p_idx);}
   void setUid(const size_t p_idx, const walberla::id_t& v) { ps_->setUid(p_idx, v);}
   
   const walberla::id_t& getId1(const size_t p_idx) const {return ps_->getId1(p_idx);}
   walberla::id_t& getId1Ref(const size_t p_idx) {return ps_->getId1Ref(p_idx);}
   void setId1(const size_t p_idx, const walberla::id_t& v) { ps_->setId1(p_idx, v);}
   
   const walberla::id_t& getId2(const size_t p_idx) const {return ps_->getId2(p_idx);}
   walberla::id_t& getId2Ref(const size_t p_idx) {return ps_->getId2Ref(p_idx);}
   void setId2(const size_t p_idx, const walberla::id_t& v) { ps_->setId2(p_idx, v);}
   
   const real_t& getDistance(const size_t p_idx) const {return ps_->getDistance(p_idx);}
   real_t& getDistanceRef(const size_t p_idx) {return ps_->getDistanceRef(p_idx);}
   void setDistance(const size_t p_idx, const real_t& v) { ps_->setDistance(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getNormal(const size_t p_idx) const {return ps_->getNormal(p_idx);}
   walberla::mesa_pd::Vec3& getNormalRef(const size_t p_idx) {return ps_->getNormalRef(p_idx);}
   void setNormal(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setNormal(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const {return ps_->getPosition(p_idx);}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t p_idx) {return ps_->getPositionRef(p_idx);}
   void setPosition(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setPosition(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getT(const size_t p_idx) const {return ps_->getT(p_idx);}
   walberla::mesa_pd::Vec3& getTRef(const size_t p_idx) {return ps_->getTRef(p_idx);}
   void setT(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setT(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getO(const size_t p_idx) const {return ps_->getO(p_idx);}
   walberla::mesa_pd::Vec3& getORef(const size_t p_idx) {return ps_->getORef(p_idx);}
   void setO(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setO(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getR1(const size_t p_idx) const {return ps_->getR1(p_idx);}
   walberla::mesa_pd::Vec3& getR1Ref(const size_t p_idx) {return ps_->getR1Ref(p_idx);}
   void setR1(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setR1(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getR2(const size_t p_idx) const {return ps_->getR2(p_idx);}
   walberla::mesa_pd::Vec3& getR2Ref(const size_t p_idx) {return ps_->getR2Ref(p_idx);}
   void setR2(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setR2(p_idx, v);}
   
   const real_t& getMu(const size_t p_idx) const {return ps_->getMu(p_idx);}
   real_t& getMuRef(const size_t p_idx) {return ps_->getMuRef(p_idx);}
   void setMu(const size_t p_idx, const real_t& v) { ps_->setMu(p_idx, v);}
   
   const walberla::mesa_pd::Vec3& getP(const size_t p_idx) const {return ps_->getP(p_idx);}
   walberla::mesa_pd::Vec3& getPRef(const size_t p_idx) {return ps_->getPRef(p_idx);}
   void setP(const size_t p_idx, const walberla::mesa_pd::Vec3& v) { ps_->setP(p_idx, v);}
   
   const walberla::mesa_pd::Mat3& getDiag_nto(const size_t p_idx) const {return ps_->getDiag_nto(p_idx);}
   walberla::mesa_pd::Mat3& getDiag_ntoRef(const size_t p_idx) {return ps_->getDiag_ntoRef(p_idx);}
   void setDiag_nto(const size_t p_idx, const walberla::mesa_pd::Mat3& v) { ps_->setDiag_nto(p_idx, v);}
   
   const walberla::mesa_pd::Mat3& getDiag_nto_inv(const size_t p_idx) const {return ps_->getDiag_nto_inv(p_idx);}
   walberla::mesa_pd::Mat3& getDiag_nto_invRef(const size_t p_idx) {return ps_->getDiag_nto_invRef(p_idx);}
   void setDiag_nto_inv(const size_t p_idx, const walberla::mesa_pd::Mat3& v) { ps_->setDiag_nto_inv(p_idx, v);}
   
   const walberla::mesa_pd::Mat2& getDiag_to_inv(const size_t p_idx) const {return ps_->getDiag_to_inv(p_idx);}
   walberla::mesa_pd::Mat2& getDiag_to_invRef(const size_t p_idx) {return ps_->getDiag_to_invRef(p_idx);}
   void setDiag_to_inv(const size_t p_idx, const walberla::mesa_pd::Mat2& v) { ps_->setDiag_to_inv(p_idx, v);}
   
   const real_t& getDiag_n_inv(const size_t p_idx) const {return ps_->getDiag_n_inv(p_idx);}
   real_t& getDiag_n_invRef(const size_t p_idx) {return ps_->getDiag_n_invRef(p_idx);}
   void setDiag_n_inv(const size_t p_idx, const real_t& v) { ps_->setDiag_n_inv(p_idx, v);}
   

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
   const walberla::id_t& getUid(const size_t /*p_idx*/) const {return uid_;}
   void setUid(const size_t /*p_idx*/, const walberla::id_t& v) { uid_ = v;}
   walberla::id_t& getUidRef(const size_t /*p_idx*/) {return uid_;}
   
   const walberla::id_t& getId1(const size_t /*p_idx*/) const {return id1_;}
   void setId1(const size_t /*p_idx*/, const walberla::id_t& v) { id1_ = v;}
   walberla::id_t& getId1Ref(const size_t /*p_idx*/) {return id1_;}
   
   const walberla::id_t& getId2(const size_t /*p_idx*/) const {return id2_;}
   void setId2(const size_t /*p_idx*/, const walberla::id_t& v) { id2_ = v;}
   walberla::id_t& getId2Ref(const size_t /*p_idx*/) {return id2_;}
   
   const real_t& getDistance(const size_t /*p_idx*/) const {return distance_;}
   void setDistance(const size_t /*p_idx*/, const real_t& v) { distance_ = v;}
   real_t& getDistanceRef(const size_t /*p_idx*/) {return distance_;}
   
   const walberla::mesa_pd::Vec3& getNormal(const size_t /*p_idx*/) const {return normal_;}
   void setNormal(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { normal_ = v;}
   walberla::mesa_pd::Vec3& getNormalRef(const size_t /*p_idx*/) {return normal_;}
   
   const walberla::mesa_pd::Vec3& getPosition(const size_t /*p_idx*/) const {return position_;}
   void setPosition(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { position_ = v;}
   walberla::mesa_pd::Vec3& getPositionRef(const size_t /*p_idx*/) {return position_;}
   
   const walberla::mesa_pd::Vec3& getT(const size_t /*p_idx*/) const {return t_;}
   void setT(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { t_ = v;}
   walberla::mesa_pd::Vec3& getTRef(const size_t /*p_idx*/) {return t_;}
   
   const walberla::mesa_pd::Vec3& getO(const size_t /*p_idx*/) const {return o_;}
   void setO(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { o_ = v;}
   walberla::mesa_pd::Vec3& getORef(const size_t /*p_idx*/) {return o_;}
   
   const walberla::mesa_pd::Vec3& getR1(const size_t /*p_idx*/) const {return r1_;}
   void setR1(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { r1_ = v;}
   walberla::mesa_pd::Vec3& getR1Ref(const size_t /*p_idx*/) {return r1_;}
   
   const walberla::mesa_pd::Vec3& getR2(const size_t /*p_idx*/) const {return r2_;}
   void setR2(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { r2_ = v;}
   walberla::mesa_pd::Vec3& getR2Ref(const size_t /*p_idx*/) {return r2_;}
   
   const real_t& getMu(const size_t /*p_idx*/) const {return mu_;}
   void setMu(const size_t /*p_idx*/, const real_t& v) { mu_ = v;}
   real_t& getMuRef(const size_t /*p_idx*/) {return mu_;}
   
   const walberla::mesa_pd::Vec3& getP(const size_t /*p_idx*/) const {return p_;}
   void setP(const size_t /*p_idx*/, const walberla::mesa_pd::Vec3& v) { p_ = v;}
   walberla::mesa_pd::Vec3& getPRef(const size_t /*p_idx*/) {return p_;}
   
   const walberla::mesa_pd::Mat3& getDiag_nto(const size_t /*p_idx*/) const {return diag_nto_;}
   void setDiag_nto(const size_t /*p_idx*/, const walberla::mesa_pd::Mat3& v) { diag_nto_ = v;}
   walberla::mesa_pd::Mat3& getDiag_ntoRef(const size_t /*p_idx*/) {return diag_nto_;}
   
   const walberla::mesa_pd::Mat3& getDiag_nto_inv(const size_t /*p_idx*/) const {return diag_nto_inv_;}
   void setDiag_nto_inv(const size_t /*p_idx*/, const walberla::mesa_pd::Mat3& v) { diag_nto_inv_ = v;}
   walberla::mesa_pd::Mat3& getDiag_nto_invRef(const size_t /*p_idx*/) {return diag_nto_inv_;}
   
   const walberla::mesa_pd::Mat2& getDiag_to_inv(const size_t /*p_idx*/) const {return diag_to_inv_;}
   void setDiag_to_inv(const size_t /*p_idx*/, const walberla::mesa_pd::Mat2& v) { diag_to_inv_ = v;}
   walberla::mesa_pd::Mat2& getDiag_to_invRef(const size_t /*p_idx*/) {return diag_to_inv_;}
   
   const real_t& getDiag_n_inv(const size_t /*p_idx*/) const {return diag_n_inv_;}
   void setDiag_n_inv(const size_t /*p_idx*/, const real_t& v) { diag_n_inv_ = v;}
   real_t& getDiag_n_invRef(const size_t /*p_idx*/) {return diag_n_inv_;}
   

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
   walberla::id_t uid_;
   walberla::id_t id1_;
   walberla::id_t id2_;
   real_t distance_;
   walberla::mesa_pd::Vec3 normal_;
   walberla::mesa_pd::Vec3 position_;
   walberla::mesa_pd::Vec3 t_;
   walberla::mesa_pd::Vec3 o_;
   walberla::mesa_pd::Vec3 r1_;
   walberla::mesa_pd::Vec3 r2_;
   real_t mu_;
   walberla::mesa_pd::Vec3 p_;
   walberla::mesa_pd::Mat3 diag_nto_;
   walberla::mesa_pd::Mat3 diag_nto_inv_;
   walberla::mesa_pd::Mat2 diag_to_inv_;
   real_t diag_n_inv_;
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla