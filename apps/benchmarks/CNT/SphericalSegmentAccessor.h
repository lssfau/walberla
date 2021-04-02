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
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/ParticleAccessor.h"
#include "mesa_pd/kernel/cnt/Parameters.h"

namespace walberla {
namespace mesa_pd {

class SphericalSegmentAccessor : public data::ParticleAccessor
{
public:
   SphericalSegmentAccessor(std::shared_ptr<data::ParticleStorage>& ps)
         : ParticleAccessor(ps)
   {}

   constexpr auto getMass(const size_t /*p_idx*/) const {return kernel::cnt::mass_T;}
   walberla::real_t& getMassRef(const size_t /*p_idx*/) {WALBERLA_ABORT("not available");}
   void setMass(const size_t /*p_idx*/, walberla::real_t const & /*v*/) {WALBERLA_ABORT("not available");}

   constexpr auto getInvMass(const size_t /*p_idx*/) const {return 1_r / kernel::cnt::mass_T;}
   walberla::real_t& getInvMassRef(const size_t /*p_idx*/) {WALBERLA_ABORT("not available");}
   void setInvMass(const size_t /*p_idx*/, walberla::real_t const & /*v*/) {WALBERLA_ABORT("not available");}

   constexpr auto& getInertiaBF(const size_t /*p_idx*/) const {return Ia;}
   Mat3& getInertiaBFRef(const size_t /*p_idx*/) {WALBERLA_ABORT("not available");}
   void setInertiaBF(const size_t /*p_idx*/, Mat3 const & /*v*/) {WALBERLA_ABORT("not available");}

   constexpr auto& getInvInertiaBF(const size_t /*p_idx*/) const {return invI;}
   Mat3& getInvInertiaBFRef(const size_t /*p_idx*/) {WALBERLA_ABORT("not available");}
   void setInvInertiaBF(const size_t /*p_idx*/, Mat3 const & /*v*/) {WALBERLA_ABORT("not available");}
private:
   //  - sphere   :  I = (2/5)*mass*radius^2
   static constexpr auto Ia = 0.4_r * kernel::cnt::mass_T * kernel::cnt::inner_radius * kernel::cnt::inner_radius;
   static constexpr auto invI = Mat3(1_r/Ia, 0_r, 0_r, 0_r, 1_r/Ia, 0_r, 0_r, 0_r, 1_r/Ia);
};

} //namespace mesa_pd
} //namespace walberla
