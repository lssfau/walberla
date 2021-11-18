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
//! \file ContactResolvers.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for DEM contact models
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "ICR.h"
#include "pe/contact/Contact.h"
#include "pe/contact/ContactFunctions.h"

namespace walberla {
namespace pe {
namespace cr {

class ResolveContactSpringDashpotHaffWerner {
public:
   void operator()( ContactID c, const real_t /*dt*/ ) const
   {
      WALBERLA_LOG_DETAIL( "resolving contact: " << c->getID() );
      BodyID b1( c->getBody1()->getTopSuperBody() );
      BodyID b2( c->getBody2()->getTopSuperBody() );

      // Global position of contact
      const Vec3 gpos( c->getPosition() );

      // The absolute value of the penetration length
      real_t delta( -c->getDistance() );

      // Calculating the relative velocity in normal and tangential direction
      // The negative signs result from the different definition of the relative
      // normal velocity of the pe (see Contact::getType)
      const real_t relVelN( -c->getNormalRelVel() );
      const Vec3   relVel ( -c->getRelVel() );
      const Vec3   relVelT( relVel - ( relVelN * c->getNormal() ) );

      real_t fNabs( 0 );
      Vec3   fN;

      // Calculating the normal force based on a linear spring-dashpot force model
      fNabs = getStiffness(c) * delta + getDampingN(c) * relVelN;
      if( fNabs < real_c(0) ) fNabs = real_c(0);
      fN = fNabs * c->getNormal();

      // Calculating the tangential force based on the model by Haff and Werner
      const real_t fTabs( std::min( getDampingT(c) * relVelT.length(), getFriction(c) * fNabs ) );
      const Vec3   fT   ( fTabs * relVelT.getNormalizedIfNotZero() );

      // Add normal force at contact point
      b1->addForceAtPos(  fN, gpos );
      b2->addForceAtPos( -fN, gpos );

      // Add tangential force at contact point
      b1->addForceAtPos(  fT, gpos );
      b2->addForceAtPos( -fT, gpos );

      WALBERLA_LOG_DETAIL("nForce: " << fN << "\ntForce: " << fT);
   }
};

}  // namespace cr
} // namespace pe
}  // namespace walberla
