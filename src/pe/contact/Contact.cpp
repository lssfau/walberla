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
//! \file Contact.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the Contact class.
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include "pe/Materials.h"
#include "pe/contact/Contact.h"
#include <core/debug/Debug.h>
#include "core/logging/all.h"


namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Contact constructor for a vertex/face contact.
 *
 * \param g1 The first contacting geometric primitive.
 * \param g2 The second contacting geometric primitive.
 * \param gpos The global position of the contact.
 * \param normal The global normal of the contact (running from body 2 to body 1).
 * \param dist The distance between the surfaces of the contacting rigid bodies.
 */
Contact::Contact( GeomID g1, GeomID g2, const Vec3& gpos, const Vec3& normal, real_t dist )
   : b1_( g1 )  // The first contacting rigid body
   , b2_( g2 )  // The second contacting rigid body
   , gpos_(gpos)                   // Global position
   , normal_(normal)               // Normal of the contact
   , e1_()                         // Edge direction of the colliding edge of body 1
   , e2_()                         // Edge direction of the colliding edge of body 2
   , dist_(dist)                   // Distance between the surfaces of the contacting rigid bodies
{
   WALBERLA_ASSERT_FLOAT_EQUAL( normal.sqrLength(), real_c(1), "Invalid contact normal\n" << g1 << "\n" << g2 );
   WALBERLA_ASSERT( !( b1_->hasInfiniteMass() && b2_->hasInfiniteMass() ), "Invalid contact between two rigid bodies with infinite masses" );

   // Debugging output
   WALBERLA_LOG_DETAIL( "         => contact-id = " << id_ );
}
//*************************************************************************************************

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the normal relative acceleration between the two contacting rigid bodies.
 *
 * \return The relative acceleration in normal direction.
 */
//real_t Contact::getNormalRelAcc() const
//{
//   return normal_ * ( b1_->accFromWF( gpos_ ) - b2_->accFromWF( gpos_ ) ) +
//          real_c(2) *  getNDot()  * ( b1_->velFromWF( gpos_ ) - b2_->velFromWF( gpos_ ) );
//}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Global output operator for ContactType.
 *
 * \param os Reference to the output stream.
 * \param type The ContactType to be put into the stream.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, Contact::ContactType type )
{
   switch( type )
   {
      case Contact::colliding:  os << "colliding";  break;
      case Contact::resting:    os << "resting";    break;
      case Contact::separating: os << "separating"; break;
      default: WALBERLA_ASSERT( false, "Unknown contact type" ); break;
   }

   return os;
}
//*************************************************************************************************

} // namespace pe
}  // namespace walberla
