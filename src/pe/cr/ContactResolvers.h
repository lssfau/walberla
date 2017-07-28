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
   void operator()( ContactID c ) const
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
      const Vec3   fT   ( fTabs * relVelT.getNormalizedOrZero() );

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
