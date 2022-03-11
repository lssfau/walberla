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
//! \author Christoph Rettinger <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class ParticleAccessorWithBaseShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithBaseShape(std::shared_ptr<data::ParticleStorage>& ps)
      : ParticleAccessor(ps)
   {}

   const auto& getInvMass(const size_t p_idx) const {return ps_->getBaseShapeRef(p_idx)->getInvMass();}
   const auto& getMass(const size_t p_idx) const {return ps_->getBaseShapeRef(p_idx)->getMass();}

   const auto& getInvInertiaBF(const size_t p_idx) const {return ps_->getBaseShapeRef(p_idx)->getInvInertiaBF();}
   const auto& getInertiaBF(const size_t p_idx) const {return ps_->getBaseShapeRef(p_idx)->getInertiaBF();}

   auto getVolume(const size_t p_idx) const {return ps_->getBaseShapeRef(p_idx)->getVolume();}

   data::BaseShape* getShape(const size_t p_idx) const {return ps_->getBaseShape(p_idx).get();}
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla
