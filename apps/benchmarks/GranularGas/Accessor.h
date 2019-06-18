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
//! \file   Accessor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

namespace walberla {
namespace mesa_pd {

class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss)
      : ParticleAccessor(ps)
      , ss_(ss)
   {}

   const walberla::real_t& getInvMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass();}
   walberla::real_t& getInvMassRef(const size_t p_idx) {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass();}
   void setInvMass(const size_t p_idx, const walberla::real_t& v) { ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass() = v;}

   const auto& getInvInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF();}
   auto& getInvInertiaBFRef(const size_t p_idx) {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF();}
   void setInvInertiaBF(const size_t p_idx, const Mat3& v) { ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF() = v;}

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};

} // namespace mesa_pd
} // namespace walberla
