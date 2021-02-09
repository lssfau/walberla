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
//! \file ICR.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/BodyStorage.h"
#include "core/DataTypes.h"
#include "core/timing/TimingTree.h"

namespace walberla {
namespace pe {
namespace cr {

class ICR {
public:
   ICR() : /*globalLinearDrag_(0),*/ globalLinearAcceleration_(0) {}
   virtual ~ICR() = default;

   virtual void timestep( const real_t dt ) = 0;

//   /// Sets the global drag coefficient.
//   /// During the integration of the linear velocity of a rigid body
//   /// a force equal to \code - drag * RigidBody::getLinearVel() \endcode is applied!
//   /// \attention In contrast to the former Settings::damping() this will result in different
//   /// final positions when the timestep length is changed!
//   inline void        setGlobalLinearDrag(const real_t drag) { globalLinearDrag_ = drag;}
//   /// Gets the global drag coefficient.
//   /// \see setGlobalLinearDrag
//   inline real_t      getGlobalLinearDrag() { return globalLinearDrag_; }

   /// Sets the global linear acceleration.
   /// This can be used for example to set a gravitational force.
   inline void        setGlobalLinearAcceleration(const Vec3& acc) { globalLinearAcceleration_ = acc; }
   inline const Vec3& getGlobalLinearAcceleration() const { return globalLinearAcceleration_; }

   virtual inline real_t            getMaximumPenetration()        const;
   virtual inline size_t            getNumberOfContacts()          const;
   virtual inline size_t            getNumberOfContactsTreated()   const;
private:
//   real_t globalLinearDrag_;
   Vec3   globalLinearAcceleration_;
};

inline real_t ICR::getMaximumPenetration()        const
{
   static bool warned = false;
   if (!warned) {
      WALBERLA_LOG_WARNING("getMaximumPenetration() is not implemented for this solver!");
      warned = true;
   }

   return real_c(0);
}
inline size_t ICR::getNumberOfContacts()          const
{
   static bool warned = false;
   if (!warned) {
      WALBERLA_LOG_WARNING("getMaximumPenetration() is not implemented for this solver!");
      warned = true;
   }

   return 0;
}
inline size_t ICR::getNumberOfContactsTreated()   const
{
   static bool warned = false;
   if (!warned) {
      WALBERLA_LOG_WARNING("getMaximumPenetration() is not implemented for this solver!");
      warned = true;
   }

   return 0;
}

}  // namespace cr
}  // namespace pe
}  // namespace walberla
