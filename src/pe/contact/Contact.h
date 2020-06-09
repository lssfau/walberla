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
//! \file Contact.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the Contact class.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/GeomPrimitive.h>
#include <pe/Thresholds.h>
#include <pe/Types.h>
#include <core/math/Vector3.h>
#include <pe/Config.h>
#include <core/NonCopyable.h>
#include <core/DataTypes.h>


namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Contact between rigid bodies.
 *
 * The Contact class is the base class for all types of contacts between rigid bodies in the
 * simulation system. Contacts between rigid bodies are classified depending on the relative
 * velocity between the touching rigid bodies:
 *
 *  - \f$ v_{rel} \f$ > 0: separating contact
 *  - \f$ v_{rel} \f$ = 0: resting contact
 *  - \f$ v_{rel} \f$ < 0: colliding contact\n
 *
 * (in order to classify the contact, pe::collisionThreshold is used as tolerance level).
 */
class Contact
{
public:
   //=================================================================================================
   //
   //  CONTACT TYPES
   //
   //=================================================================================================

   //*************************************************************************************************
   //! Classification of contacts.
   /*! Contacts between rigid bodies are classified depending on the relative velocity between
    *  the touching rigid bodies:
    *
    *   - \f$ v_{rel} \f$ > 0: separating contact
    *   - \f$ v_{rel} \f$ = 0: resting contact
    *   - \f$ v_{rel} \f$ < 0: colliding contact\n
    *
    * (in order to classify the contact, pe::collisionThreshold is used as tolerance level).
    */
   enum ContactType {
      colliding  = 0,  //!< Colliding contacts (vrel < 0).
      resting    = 1,  //!< Resting contacts (vrel = 0).
      separating = 2   //!< Separating contacts (vrel > 0).
   };
   //*************************************************************************************************
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   Contact( GeomID g1, GeomID g2, const Vec3& gpos, const Vec3& normal, real_t dist );
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline bool         isPenetrating()                  const;
   inline bool         isContacting( ConstBodyID body ) const;

   inline id_t           getID()             const;
   inline size_t         getIndex()          const;
   inline GeomID         getBody1()          const;
   inline GeomID         getBody2()          const;
   inline const Vec3&    getPosition()       const;
   inline const Vec3&    getNormal()         const;
//   inline const Vec3     getNDot()           const;
   inline real_t         getDistance()       const;

   inline ContactType    getType()           const;

   inline real_t         getNormalRelVel()   const;
   inline const Vec3     getRelVel()         const;
//          real_t         getNormalRelAcc()   const;
          real_t         getNormalRelForce() const;
   //@}
   //**********************************************************************************************

protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   id_t id_;           //!< User-specific contact ID.

   GeomID b1_;         //!< The first contacting rigid body.
   GeomID b2_;         //!< The second contacting rigid body.
   Vec3 gpos_;         //!< The global position of the contact.
   Vec3 normal_;       //!< Normal of the contact.
                       /*!< The normal is defined within the global world frame and points
                            from body 2 to body 1. */
   Vec3 e1_;           //!< Edge direction of the colliding edge of body 1.
   Vec3 e2_;           //!< Edge direction of the colliding edge of body 2.
   real_t dist_;         //!< Distance between the surfaces of the contacting rigid bodies.
                       /*!< A negative distance means penetration of the two bodies. */
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Checks if the two contacting rigid bodies are penetrating each other.
 *
 * \return \a true if the two rigid bodies are penetrating each other, \a false if not.
 *
 * The tolerance level of the check is pe::contactThreshold.
 */
inline bool Contact::isPenetrating() const
{
   return dist_ < -contactThreshold;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given rigid body is involved in this contact.
 *
 * \return \a true if the rigid body is involved in this contact, \a false if not.
 */
inline bool Contact::isContacting( ConstBodyID body ) const
{
   return ( b1_ == body || b2_ == body );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the ID of the contact.
 *
 * \return The contact ID.
 */
inline id_t Contact::getID() const
{
   return id_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the first constrained rigid body.
 *
 *\return The first constrained rigid body.
 */
inline GeomID Contact::getBody1() const
{
   return b1_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the second constrained rigid body.
 *
 *\return The second constrained rigid body.
 */
inline GeomID Contact::getBody2() const
{
   return b2_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the global position of the contact.
 *
 * \return Global position of the contact.
 */
inline const Vec3& Contact::getPosition() const
{
   return gpos_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the global normal of the contact.
 *
 * \return Global normal of the contact.
 */
inline const Vec3& Contact::getNormal() const
{
   return normal_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the derivative of the normal vector.
 *
 * \return Derivative of the normal vector.
 */
//inline const Vec3 Contact::getNDot() const
//{
//   if( vf_ ) {
//      return b2_->getAngularVel() % normal_;
//   }
//   else {
//      const Vec3 e1dot( b1_->getAngularVel() % e1_ );
//      const Vec3 e2dot( b2_->getAngularVel() % e2_ );
//      const Vec3 z( e1dot % e2_ + e1_ % e2dot );

//      return z - ( ( z % normal_ ) % normal_ );
//   }
//}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the distance between the surfaces of the contacting rigid bodies.
 *
 * \return Distance between the surfaces of the contacting rigid bodies.
 */
inline real_t Contact::getDistance() const
{
   return dist_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the type of the contact.
 *
 * \return Type of the contact.
 *
 * This function returns the type of the contact: separating, resting or colliding. The contact
 * is classified according to the relative velocity of the two touching rigid bodies in the
 * contact:
 *
 *  - \f$ v_{rel} \f$ > 0: separating contact
 *  - \f$ v_{rel} \f$ = 0: resting contact
 *  - \f$ v_{rel} \f$ < 0: colliding contact\n
 *
 * (in order to classify the contact, pe::collisionThreshold is used as tolerance level).
 */
inline Contact::ContactType Contact::getType() const
{
   real_t vrel = normal_ * ( b1_->velFromWF( gpos_ ) - b2_->velFromWF( gpos_ ) );

   if( vrel > collisionThreshold ) return Contact::separating;
   else if( vrel > -collisionThreshold ) return Contact::resting;
   else return Contact::colliding;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the normal relative velocity between the two contacting rigid bodies.
 *
 * \return The relative velocity in normal direction.
 */
inline real_t Contact::getNormalRelVel() const
{
   return normal_ * ( b1_->velFromWF( gpos_ ) - b2_->velFromWF( gpos_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the relative velocity between the two contacting rigid bodies.
 *
 * \return The relative velocity.
 */
inline const Vec3 Contact::getRelVel() const
{
   return b1_->velFromWF( gpos_ ) - b2_->velFromWF( gpos_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for contacts.
 *
 * \param os Reference to the output stream.
 * \param c Reference to a constant contact object.
 * \return Reference to the output stream.
 */
inline std::ostream& operator<<( std::ostream& os, const Contact& c )
{
   const char* tab="";
   os << "--" << "CONTACT PARAMETERS"
      << "---------------------------------------------------------\n";
   using std::setw;

   os << tab << " Contact " << c.getID() << " at position " << c.getPosition() << "\n";

   os << tab << "   Body 1            = " << c.getBody1()->getSystemID() << "\n"
      << tab << "   Body 2            = " << c.getBody2()->getSystemID() << "\n"
      << tab << "   Normal            = " << c.getNormal() << "\n"
      << tab << "   Distance          = " << c.getDistance() << "\n"
      << tab << "   Locality          = " << ( c.getBody1()->isRemote() ? "remote" : "local") << "-" << ( c.getBody2()->isRemote() ? "remote" : "local") << " contact\n";

   os << tab << "   Contact type      = " << (c.getType() == Contact::colliding ? "colliding" : (c.getType() == Contact::resting ? "resting" : "separating")) << "\n"
      << tab << "   Rel. normal vel.  = " << c.getNormalRelVel() << "\n";
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for contact handles.
 *
 * \param os Reference to the output stream.
 * \param c Constant contact handle.
 * \return Reference to the output stream.
 */
inline std::ostream& operator<<( std::ostream& os, ConstContactID c )
{
   return operator<<( os, *c );
}
//*************************************************************************************************

//=================================================================================================
//
//  CONTACT TYPE UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Contact type utility functions */
//@{
std::ostream& operator<<( std::ostream& os, Contact::ContactType type );
//@}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
