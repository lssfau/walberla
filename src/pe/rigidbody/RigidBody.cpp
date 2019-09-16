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
//! \file RigidBody.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "RigidBody.h"

namespace walberla{
namespace pe{

//*************************************************************************************************
/*!\brief Constructor for the RigidBody class.
 *
 * \param finite Specifies if the rigid body is finite or not.
 * \param visible Specifies if the rigid body is visible or not.
 * \param sid The unique system-specific ID of the rigid body.
 * \param uid The user-specific ID of the rigid body.
 */
RigidBody::RigidBody( id_t const typeID, id_t sid, id_t uid )
   : awake_( true )           // Sleep mode flag
   , mass_( 0 )               // Total mass of the rigid body
   , invMass_( 0 )            // Inverse total mass of the rigid body
   , motion_(sleepThreshold)  // The current motion of the rigid body
   , v_()                     // Linear velocity
   , w_()                     // Angular velocity
   , force_()                 // Total force
   , torque_()                // Total torque
   , I_( real_c(0) )          // Moment of inertia
   , Iinv_( real_c(0) )       // Inverse moment of inertia
   , manager_(nullptr)        // The rigid body manager responsible for the rigid body
   , finite_ (true)           // Finiteness flag
   , visible_(true)           // Visibility flag
   , remote_ (false)          // Remote flag
   , communicating_  (false)  // Local flag
   , global_ (false)          // Global flag
   , toBeDeleted_(false)      // deletion flag
   , sid_    (sid)            // System-specific body index
   , uid_    (uid)            // User-specific body ID
   , gpos_()                  // Global position of the center of mass
   , q_()                     // Orientation of the body frame
   , R_()                     // Rigid body rotation
   , typeID_(typeID)          // geometry type
{
   sb_ = this;           // The superordinate rigid body
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the RigidBody class.
 */
RigidBody::~RigidBody()
= default;
//*************************************************************************************************




//=================================================================================================
//
//  DEBUGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks the validity of the state of the rigid body.
 *
 * \return \a true if no error is detected, \a false if any error if found.
 */
bool RigidBody::checkInvariants()
{
   bool error( false );

   // Checking that an infinite rigid body is fixed
   if( !isFinite() && !isFixed() ) {
      std::cerr << "\n **** Infinite rigid body is not fixed ****";
      error = true;
   }

   // Checking that a global rigid body is local
   if( isGlobal() && isRemote() ) {
      std::cerr << "\n **** Global rigid body is not local ****";
      error = true;
   }

   // Checking the mass properties
   if( mass_ < real_c(0) || invMass_ < real_c(0) ) {
      std::cerr << "\n **** Invalid mass properties detected ****";
      error = true;
   }

   // Checking the mass properties of a fixed rigid body
   if( isFixed() && ( invMass_ > real_c(0) || Iinv_ != Mat3(real_c(0)) ) ) {
      std::cerr << "\n **** Invalid mass properties for fixed rigid body ****";
      error = true;
   }

   // Checking that a fixed or sleeping rigid body has zero linear and angular velocity
#if MOBILE_INFINITE
   if( !isAwake() && ( v_ != real_c(0) || w_ != real_c(0) ) ) {
#else
   if( ( isFixed() || !isAwake() ) && ( v_ != real_c(0) || w_ != real_c(0) ) ) {
#endif
      std::cerr << "\n **** Invalid velocity for immobile rigid body " << uid_ << " ****";
      error = true;
   }

   // Checking for nan-values
   if( isnan(  gpos_  ) || isnan( v_ ) || isnan( w_ ) || isnan( force_ ) ||
       isnan( torque_ ) || isnan( q_ ) || isnan( R_ ) ) {
      std::cerr << "\n **** Nan-value detected in rigid body ****";
      error = true;
   }

   // Printing the current state of the rigid body
   if( error ) {
      std::cerr << std::boolalpha << "\n"
                << "   User-ID           = " << getID() << "\n"
                << "   System-ID         = " << getSystemID() << "\n"
                << "   Finite            = " << isFinite() << "\n"
                << "   Awake             = " << isAwake() << "\n"
                << "   Fixed             = " << isFixed() << "\n"
                << "   Remote            = " << isRemote() << "\n"
                << "   Communicating     = " << isCommunicating() << "\n"
                << "   Total mass        = " << getMass() << "\n"
                << "   Inverse mass      = " << getInvMass() << "\n"
                << "   Global position   = " << getPosition() << "\n"
                << "   Relative position = " << getRelPosition() << "\n"
                << "   Linear velocity   = " << getLinearVel() << "\n"
                << "   Angular velocity  = " << getAngularVel() << "\n"
                << "   Acting force      = " << getForce() << "\n"
                << "   Acting torque     = " << getTorque() << "\n"
                << "   Bounding box      = " << getAABB() << "\n"
                << "   Quaternion        = " << getQuaternion() << "\n"
                << "   Rotation matrix   = ( " << R_[0] << " , " << R_[1] << " , " << R_[2] << " )\n"
                << "                       ( " << R_[3] << " , " << R_[4] << " , " << R_[5] << " )\n"
                << "                       ( " << R_[6] << " , " << R_[7] << " , " << R_[8] << " )\n"
                << std::setiosflags(std::ios::right)
                << "   Moment of inertia :\n"
                << std::setw(18) << I_[0] << std::setw(18) << I_[1] << std::setw(18) << I_[2] << "\n"
                << std::setw(18) << I_[3] << std::setw(18) << I_[4] << std::setw(18) << I_[5] << "\n"
                << std::setw(18) << I_[6] << std::setw(18) << I_[7] << std::setw(18) << I_[8] << "\n"
                << "   Inverse moment of inertia :\n"
                << std::setw(18) << Iinv_[0] << std::setw(18) << Iinv_[1] << std::setw(18) << Iinv_[2] << "\n"
                << std::setw(18) << Iinv_[3] << std::setw(18) << Iinv_[4] << std::setw(18) << Iinv_[5] << "\n"
                << std::setw(18) << Iinv_[6] << std::setw(18) << Iinv_[7] << std::setw(18) << Iinv_[8] << "\n"
                << std::resetiosflags(std::ios::right)
                << std::endl;

      return false;
   }
   else return true;
}
//*************************************************************************************************



//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a relative point.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point relative to the body's center of mass.
 */
Vec3 RigidBody::velFromBF( real_t px, real_t py, real_t pz ) const
{
   return velFromBF( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a relative point.
 *
 * \param rpos The relative coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point relative to the body's center of mass.
 */
Vec3 RigidBody::velFromBF( const Vec3& rpos ) const
{
   if( !hasSuperBody() )
      return v_ + w_ % ( R_ * rpos );
   else return sb_->velFromWF( gpos_ + R_ * rpos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point in global coordinates.
 */
Vec3 RigidBody::velFromWF( real_t px, real_t py, real_t pz ) const
{
   return velFromWF( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point in global coordinates.
 */
Vec3 RigidBody::velFromWF( const Vec3& gpos ) const
{
   if( !hasSuperBody() )
      return v_ + w_ % ( gpos - gpos_ );
   else return sb_->velFromWF( gpos );
}
//*************************************************************************************************



//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for rigid bodies.
 *
 * \param os Reference to the output stream.
 * \param b Reference to a constant rigid body object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const RigidBody& b )
{
   os << "--" << "RIGID BODY PARAMETERS"
      << "---------------------------------------------------------\n";
   b.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for rigid body handles.
 *
 * \param os Reference to the output stream.
 * \param b Constant rigid body handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstBodyID b )
{
   os << "--" << "RIGID BODY PARAMETERS"
      << "---------------------------------------------------------\n";
   b->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
