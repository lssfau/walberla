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
//! \file ContactFunctions.impl.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "ContactFunctions.h"

#include "core/logging/Logging.h"

namespace walberla {
namespace pe {

inline
bool shouldContactBeTreated( ContactID c, const math::AABB& blkAABB )
{
   const int myRank( mpi::MPIManager::instance()->rank() );

   BodyID b1( c->getBody1() );
   BodyID b2( c->getBody2() );

   WALBERLA_ASSERT( !b1->hasInfiniteMass() || !b2->hasInfiniteMass(), "Invalid contact between two objects with infinite masses." );

   /* Contact filtering rules
    *
    * L: Local body
    * G: Global body
    * R: Remote body
    *
    * Option 1:              Option 2:
    * +---+---+---+---+      +---+---+---+---+
    * |   | L | G | R |      |   | L | G | R |
    * +---+---+---+---+      +---+---+---+---+
    * | L | + | + | * |      | L |§/+| + | § |
    * +---+---+---+---+      +---+---+---+---+
    * | G | + | ~ | - |      | G | + | ~ | - |
    * +---+---+---+---+      +---+---+---+---+
    * | R | * | - | # |      | R | § | - | § |
    * +---+---+---+---+      +---+---+---+---+
    *
    *  + Accept contact unconditionally
    *  - Reject contact unconditionally
    *  * Accept contact if we own the contact point
    *  # Accept contact if we own the contact point and the owners of the involved bodies are not the same
    *  ~ Accept contact only on root process
    *  § Accept contact if we are the owner with the smallest rank witnessing the contact or if none of the owners witness the contact and we are the process with the smallest rank witnessing the contact
    *
    * Note: Local-global contacts actually require a reduction of the contact reactions applied to the global body (unless it is fixed).
    * => MPI_Allreduce for all non-fixed global bodies before time-integration.
    */

   if( !b1->isRemote() && !b2->isRemote() ) {
      // local-local, local-global, global-global contacts

      if( b1->isGlobal() && b2->isGlobal() ) {
         // Resolve global-global contacts only on root process
         if( myRank != 0 ) {
            WALBERLA_LOG_DETAIL( "Rejecting global-global contact " << *c << " on non-root process." );
            return false;
         }
      } else
      {
         // Always resolve local-local and local-global contacts even if they are outside of our domain
      }
   }
   else {
      // local-remote, global-remote or remote-remote contact

      if( b1->isGlobal() || b2->isGlobal() ) {
         // Never resolve remote-global contacts
         WALBERLA_LOG_DETAIL( "Rejecting global-remote contact " << *c << "." );
         return false;
      }
      else if( b1->isRemote() && b2->isRemote() && b1->MPITrait.getOwner().blockID_ == b2->MPITrait.getOwner().blockID_ ) {
         WALBERLA_LOG_DETAIL( "Rejecting remote-remote contact since it will be a local-local contact at the owner process: " << *c << "." );
         return false;
      } else
      {
         // Option 1
         if( !blkAABB.contains( c->getPosition() ) )
         {
            if( b1->isRemote() && b2->isRemote() )
            {
               WALBERLA_LOG_DETAIL( "Rejecting remote-remote contact " << *c << " since we don't own it." );
               return false;
            } else
            {
               WALBERLA_LOG_DETAIL( "Rejecting remote-local contact " << *c << " since we don't own it." );
               return false;
            }
         }
      }
   }

   return true;
}

// Calculating the coefficient of restitution between the two colliding (sub-)bodies
// In case the normal relative velocity between the two colliding rigid bodies is smaller than the
// restitution threshold, a coefficient of restitution of 0 is used to prevent an infinite number
// of collisions during a single time step.
inline real_t         getRestitution(ConstContactID c)
{
   // Calculating the relative velocity
   const Vec3 rvel( c->getBody1()->velFromWF( c->getPosition() ) - c->getBody2()->velFromWF( c->getPosition() ) );  // Relative velocity
   const real_t nvel( c->getNormal() * rvel );                              // Normal relative velocity
   if( std::fabs( nvel ) > restitutionThreshold )
   {
      return Material::getRestitution( c->getBody1()->getMaterial(), c->getBody1()->getMaterial() );
   }
   return real_t(0.0);
}

// Calculate the stiffness and damping parameters
inline real_t         getStiffness(ConstContactID c)
{
   return Material::getStiffness( c->getBody1()->getMaterial(), c->getBody2()->getMaterial() );
}
inline real_t         getDampingN(ConstContactID c)
{
   return Material::getDampingN ( c->getBody1()->getMaterial(), c->getBody2()->getMaterial() );
}
inline real_t         getDampingT(ConstContactID c)
{
   return Material::getDampingT ( c->getBody1()->getMaterial(), c->getBody2()->getMaterial() );
}

// Calculating the coefficient of friction
// In case the tangential relative velocity between the two colliding rigid bodies is smaller than
// the friction threshold, the coefficient of static friction is used. Otherwise, the coefficient
// of dynamic friction is used.
inline real_t         getFriction(ConstContactID c)
{
   // Calculating the relative velocity
   const Vec3 rvel( c->getBody1()->velFromWF( c->getPosition() ) - c->getBody2()->velFromWF( c->getPosition() ) );  // Relative velocity
   const Vec3 nvel( ( c->getNormal() * rvel ) * c->getNormal() );    // Normal relative velocity
   const Vec3 tvel( rvel - nvel );                                   // Tangential relative velocity

   if( std::fabs( tvel.length() ) > frictionThreshold )
      return Material::getDynamicFriction( c->getBody1()->getMaterial(), c->getBody2()->getMaterial() );
   else
      return Material::getStaticFriction( c->getBody1()->getMaterial(), c->getBody2()->getMaterial() );
}

} // namespace pe
} // namespace walberla
