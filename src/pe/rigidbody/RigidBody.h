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
//! \file RigidBody.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "pe/Config.h"
#include "pe/ccd/HashGridsBodyTrait.h"
#include "pe/cr/HCSITSBodyTrait.h"
#include "core/math/Matrix3.h"
#include "core/math/Quaternion.h"
#include "core/math/Vector3.h"

#include "core/NonCopyable.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "pe/rigidbody/MPIRigidBodyTrait.h"

namespace walberla{
namespace pe{

template <typename BodyTypeTuple>
class Union;

/**
 * \ingroup pe
 */
class RigidBody : public ccd::HashGridsBodyTrait
                , public cr::HCSITSBodyTrait
                , private NonCopyable
{
private:
   template <typename BodyTypeTuple>
   friend class Union;

protected:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit RigidBody( id_t const typeID, id_t sid, id_t uid );
   //@}
   //**********************************************************************************************

public:
   MPIRigidBodyTrait MPITrait;

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~RigidBody() = 0;
   //@}
   //**********************************************************************************************

   //**Sleep mode functions************************************************************************
   /*!\name Sleep mode functions */
   //@{
   inline void wake();
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline bool           hasManager()        const;
   inline ManagerID      getManager();
   inline ConstManagerID getManager()        const;
   inline bool           hasSuperBody()      const;
   inline BodyID         getSuperBody();
   inline ConstBodyID    getSuperBody()      const;
   inline BodyID         getTopSuperBody();
   inline ConstBodyID    getTopSuperBody()   const;
   virtual inline bool   hasSubBodies()      const { return false; }
   inline bool           isFinite()          const;
   inline bool           isAwake()           const;
   inline bool           isFixed()           const;
   inline bool           hasInfiniteMass()   const;
   inline bool           isVisible()         const;
   inline bool           isMarkedForDeletion() const { return toBeDeleted_; }
   inline id_t           getSystemID()       const;
   inline id_t           getID()             const;
   inline const Vec3&    getRelPosition()    const;
   inline const Vec3&    getPosition()       const;
   inline const Vec3     getRelLinearVel()   const;
   inline const Vec3&    getLinearVel()      const;
   inline const Vec3     getRelAngularVel()  const;
   inline const Vec3&    getAngularVel()     const;
   inline const Quat&    getQuaternion()     const;
   inline const Mat3&    getRotation()       const;
   inline real_t         getMass()           const;
   inline real_t         getInvMass()        const;
   virtual inline real_t getVolume()         const;
   inline const Mat3&    getBodyInertia()    const;
   inline const Mat3     getInertia()        const;
   inline const Mat3&    getInvBodyInertia() const;
   inline const Mat3     getInvInertia()     const;
   inline const AABB&    getAABB()           const;
   inline real_t         getAABBSize()       const;

   inline       real_t   getKineticEnergy()     const;
   inline       real_t   getRotationalEnergy()  const;
   inline       real_t   getEnergy()            const;

   inline const Vec3     vectorFromBFtoWF( real_t vx, real_t vy, real_t vz ) const;
   inline const Vec3     vectorFromBFtoWF( const Vec3& v )             const;
   inline const Vec3     pointFromBFtoWF ( real_t px, real_t py, real_t pz ) const;
   inline const Vec3     pointFromBFtoWF ( const Vec3& rpos )          const;
   virtual const Vec3    velFromBF       ( real_t px, real_t py, real_t pz ) const;
   virtual const Vec3    velFromBF       ( const Vec3& rpos )          const;
   inline const Vec3     vectorFromWFtoBF( real_t vx, real_t vy, real_t vz ) const;
   inline const Vec3     vectorFromWFtoBF( const Vec3& v )             const;
   inline const Vec3     pointFromWFtoBF ( real_t px, real_t py, real_t pz ) const;
   inline const Vec3     pointFromWFtoBF ( const Vec3& gpos )          const;
   virtual const Vec3    velFromWF       ( real_t px, real_t py, real_t pz ) const;
   virtual const Vec3    velFromWF       ( const Vec3& gpos )          const;
   inline const Vec3     accFromWF       ( real_t px, real_t py, real_t pz ) const;
          const Vec3     accFromWF       ( const Vec3& gpos )          const;

   inline       id_t     getTypeID() const;
   //@}
   //**********************************************************************************************

   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   inline  void setSB  ( BodyID body );
   inline  void resetSB() ;

   inline  void setFinite     ( const bool finite );
   inline  void setVisible    ( bool visible );
   inline  void setPosition   ( real_t px, real_t py, real_t pz );
   inline  void setPosition   ( const Vec3& gpos );
   inline  void setOrientation( real_t r, real_t i, real_t j, real_t k );
   inline  void setOrientation( const Quat& q );
   inline  void setMassAndInertiaToInfinity();

   inline void setRelLinearVel ( real_t vx, real_t vy, real_t vz );
   inline void setRelLinearVel ( const Vec3& lvel );
   inline void setLinearVel    ( real_t vx, real_t vy, real_t vz );
   inline void setLinearVel    ( const Vec3& lvel );
   inline void setRelAngularVel( real_t ax, real_t ay, real_t az );
   inline void setRelAngularVel( const Vec3& avel );
   inline void setAngularVel   ( real_t ax, real_t ay, real_t az );
   inline void setAngularVel   ( const Vec3& avel );

   ///Marks the rigid body for deletion during the next synchronization
   /// \attention Only works on local bodies!!!<br>
   /// Has no effect on global or shadow bodies.<br>
   /// Bodies which are non communicating have to be explicitly communicated during synchronization!
   inline void markForDeletion () { toBeDeleted_ = true; }
   //@}
   //**********************************************************************************************

   //**Translation functions***********************************************************************
   /*!\name Translation functions */
   //@{
   void translate( real_t dx, real_t dy, real_t dz );
   void translate( const Vec3& dp );
   //@}
   //**********************************************************************************************

   //**Rotation functions**************************************************************************
   /*!\name Rotation functions */
   //@{
   void rotate( real_t x, real_t y, real_t z, real_t angle );
   void rotate( const Vec3& axis, real_t angle );
   void rotate( real_t xangle, real_t yangle, real_t zangle );
   void rotate( const Vec3& euler );
   void rotate( const Quat& dq );

   void rotateAroundOrigin( real_t x, real_t y, real_t z, real_t angle );
   void rotateAroundOrigin( const Vec3& axis, real_t angle );
   void rotateAroundOrigin( real_t xangle, real_t yangle, real_t zangle );
   void rotateAroundOrigin( const Vec3& euler );
   void rotateAroundOrigin( const Quat& dq );

   void rotateAroundPoint( const Vec3& point, const Vec3& axis, real_t angle );
   void rotateAroundPoint( const Vec3& point, const Vec3& euler );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   bool containsRelPoint       ( real_t px, real_t py, real_t pz ) const;
   bool containsRelPoint       ( const Vec3& rpos )                const;
   bool containsPoint          ( real_t px, real_t py, real_t pz ) const;
   bool containsPoint          ( const Vec3& gpos )                const;
   bool isSurfaceRelPoint      ( real_t px, real_t py, real_t pz ) const;
   bool isSurfaceRelPoint      ( const Vec3& rpos )                const;
   bool isSurfacePoint         ( real_t px, real_t py, real_t pz ) const;
   bool isSurfacePoint         ( const Vec3& gpos )                const;
   virtual Vec3 support                ( const Vec3& d )                   const;
   virtual Vec3 supportContactThreshold( const Vec3& d )                   const;
   //@}
   //**********************************************************************************************

   //**Force functions*****************************************************************************
   /*!\name Force functions */
   //@{
   inline bool        hasForce () const;

   inline const Vec3& getForce () const;
   inline const Vec3& getTorque() const;

   inline void        setForce           ( const Vec3& f );
   inline void        setTorque          ( const Vec3& tau );

   inline void        addRelForce        ( real_t fx, real_t fy, real_t fz );
   inline void        addRelForce        ( const Vec3& f );
   inline void        addForce           ( real_t fx, real_t fy, real_t fz );
   inline void        addForce           ( const Vec3& f );
   inline void        addRelForceAtRelPos( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz );
   inline void        addRelForceAtRelPos( const Vec3& f, const Vec3& rpos );
   inline void        addRelForceAtPos   ( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz );
   inline void        addRelForceAtPos   ( const Vec3& f, const Vec3& gpos );
   inline void        addForceAtRelPos   ( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz );
   inline void        addForceAtRelPos   ( const Vec3& f, const Vec3& rpos );
   inline void        addForceAtPos      ( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz );
   inline void        addForceAtPos      ( const Vec3& f, const Vec3& gpos );

   inline void        addTorque( real_t tx, real_t ty, real_t tz );
   inline void        addTorque( const Vec3& t );

   virtual inline void        resetForceAndTorque();
   //@}
   //**********************************************************************************************

   //**Impulse functions***************************************************************************
   /*!\name Impulse functions */
   //@{
   inline void addImpulse( real_t jx, real_t jy, real_t jz );
   inline void addImpulse( const Vec3& j );
   inline void addImpulseAtPos( real_t jx, real_t jy, real_t jz, real_t px, real_t py, real_t pz );
   inline void addImpulseAtPos( const Vec3& j, const Vec3& p );
   //@}
   //**********************************************************************************************

   //**MPI functions*******************************************************************************
   /*!\name MPI functions */
   //@{
           inline bool isRemote () const;
   virtual inline void setRemote( bool remote );
           inline bool isGlobal () const;
           inline void setGlobal( const bool global );
           inline bool isCommunicating () const;
           inline void setCommunicating( const bool communicating );
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   virtual void print( std::ostream& os, const char* tab ) const = 0;
   //@}
   //**********************************************************************************************

   //**Debugging functions*************************************************************************
   /*!\name Debugging functions */
   //@{
   virtual bool checkInvariants();
   //@}
   //**********************************************************************************************

   //**Utility functions*******************************************************************************
   /*!\name Utility functions */
   //@{
   virtual void calcBoundingBox()                                = 0;
   //@}
   //**********************************************************************************************

   //**Sleep mode functions************************************************************************
   /*!\name Sleep mode functions */
   //@{
   inline void calcMotion();
   //@}
   //**********************************************************************************************

protected:
   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   virtual void setPositionImpl       ( real_t px, real_t py, real_t pz );
   virtual void setOrientationImpl    ( real_t r, real_t i, real_t j, real_t k );
   virtual void translateImpl         ( real_t dx, real_t dy, real_t dz );
   virtual void rotateImpl            ( const Quat& dq );
   virtual void rotateAroundOriginImpl( const Quat& dq );
   virtual void rotateAroundPointImpl ( const Vec3& point, const Quat& dq );
   virtual bool containsRelPointImpl  ( real_t px, real_t py, real_t pz ) const;
   virtual bool isSurfaceRelPointImpl ( real_t px, real_t py, real_t pz ) const;

   inline  void setMassAndInertia     ( const real_t mass, const Mat3& inertia );

   inline void calcRelPosition();
   //@}
   //**********************************************************************************************

   //**Fixation functions**************************************************************************
   /*!\name Fixation functions */
   //@{
   virtual void fix  ();
   //@}
   //**********************************************************************************************

   //**Simulation functions************************************************************************
   /*!\name Simulation functions */
   //@{
   virtual void update( const Vec3& dp );  // Translation update of a subordinate rigid body
   virtual void update( const Quat& dq );  // Rotation update of a subordinate rigid body
   //@}
   //**********************************************************************************************

   //**Functions for internal changes in compound geometries***************************************
   /*!\name Functions for internal changes in compound geometries */
   //@{
   inline void signalModification();
   inline void signalTranslation();
   inline void signalRotation();
   inline void signalFixation();

   virtual void handleModification();
   virtual void handleTranslation();
   virtual void handleRotation();
   virtual void handleFixation();
   //@}
   //**********************************************************************************************

protected:

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   bool awake_;      //!< Sleep mode flag.
                  /*!< The flag value indicates if the rigid body is currently awake
                       (\a true) or in sleep mode (\a false). */
   real_t mass_;     //!< The total mass of the rigid body.
   real_t invMass_;  //!< The inverse total mass of the rigid body.
   real_t motion_;   //!< The current motion of the rigid body.
                  /*!< If this value drops under the specified sleep threshold, the
                       rigid body will be put to sleep. */
   Vec3 gpos_;       //!< The global position of the center of mass of this rigid body.
   Vec3 rpos_;       //!< The relative position within the body frame of the superordinate body.
                  /*!< If the body is contained within a superordinate Union the relative
                       position gives the position of the body's center of mass within the
                       body frame of the superordinate body. If the body is not contained
                       in a Union, the relative position is 0. */
   mutable Vec3 v_;  //!< The linear velocity of this rigid body.
   mutable Vec3 w_;  //!< Angular velocity of this rigid body.
   Vec3 force_;      //!< Total force (external+contact) acting in the body's center of mass.
   Vec3 torque_;     //!< Total torque (external+contact) acting in the body's center of mass.
   Mat3 I_;          //!< The moment of inertia in reference to the body's own body frame.
                  /*!< The moment of inertia quantifies the rotational inertia of a rigid
                       body, i.e. its inertia with respect to rotational motion, in a manner
                       somewhat analogous to how mass quantifies the inertia of a body with
                       respect to translational motion. The naming convention of the tensor
                       elements is
                                      \f[\left(\begin{array}{*{3}{c}}
                                      I_{xx} & I_{xy} & I_{xz} \\
                                      I_{yx} & I_{yy} & I_{yz} \\
                                      I_{zx} & I_{zy} & I_{zz} \\
                                      \end{array}\right)\f] */
   Mat3 Iinv_;       //!< The inverse moment of inertia within the body frame.
   Quat q_;          //!< The orientation of the body frame in the global world frame.
   Mat3 R_;          //!< The rotation in reference to the global frame of reference.
   //@}
   //**********************************************************************************************

   /*!\name Member variables */
   //@{
   ManagerID manager_;        //!< The rigid body manager responsible for the rigid body.
   BodyID sb_;                //!< The superordinate rigid body.
                              /*!< This data member is the connection to the superordinate body,
                                   which is either the enclosing Union or the rigid body itself. */
   bool finite_;              //!< Finiteness flag.
                              /*!< The flag value indicates if the rigid body is finite (\a true)
                                   or infinite (\a false). */
   bool visible_;             //!< Visibility flag.
   bool remote_;              //!< Remote flag.
                              /*!< This flag indicates whether the rigid body belongs to a remote
                                   process (\a true) or to the local process (\a false). */
   bool communicating_;       //!< Communicating flag.
                              /*!< The local flag indicates whether the rigid body is local in
                                   an MPI parallel simulation. Local bodies are not participating
                                   in any communication process. */

   bool global_;              //!< Global flag.
                              /*!< The global flag indicates whether the rigid body is global in
                                   an MPI parallel simulation. */
   bool toBeDeleted_;         //!< This flag marks the body for deletion during the next synchronization (only works on local bodies)
   id_t sid_;                 //!< The unique system-specific body ID.
   id_t uid_;                 //!< The user-specific body ID.
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   AABB aabb_;  //!< Axis-aligned bounding box for the rigid body.
   //@}
   //**********************************************************************************************
private:
   id_t typeID_; //< identify geometry type
};

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Rigid body operators */
//@{
std::ostream& operator<<( std::ostream& os, const RigidBody& b );
std::ostream& operator<<( std::ostream& os, ConstBodyID b );
//@}
//*************************************************************************************************

//=================================================================================================
//
//  SLEEP MODE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Waking the rigid body and ending the sleep mode.
 *
 * \return void
 *
 * This function wakes a rigid body from sleep mode. Note that this function has no effect if
 * it is called on a subordinate rigid body, i.e. a body contained in another rigid body.
 */
inline void RigidBody::wake()
{
   awake_ = true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating the current motion of a rigid body.
 *
 * \return void
 *
 * This function calculates the current motion of a rigid body. The motion is a scalar value,
 * consisting of the current linear and angular velocity and is calculated by the following
 * equation:

   \f[ \Phi = \mbox{bias} \cdot \Phi + (1-\mbox{bias}) \cdot ( \vec{v}^2 + \vec{w}^2 ), \f]

 * where \f$ \Phi \f$ is the motion, \f$ \vec{v} \f$ is the linear velocity, \f$ \vec{w} \f$
 * is the angular velocity and \a bias is the weighting factor for the recency-weighted average
 * between the new motion value and the motion value from the last time frame (see pe::sleepBias).
 * If the motion drops below pe::sleepThreshold, the rigid body is put to sleep.
 */
inline void RigidBody::calcMotion()
{
   WALBERLA_ASSERT( sleepThreshold >= real_c(0), "Invalid sleepThreshold value" );
   WALBERLA_ASSERT( sleepBias >= real_c(0) && sleepBias <= real_c(1), "Invalid sleepBias value" );

   // Calculating the current motion of the rigid body (recency-weighted average)
   motion_ = ( sleepBias * motion_ ) + ( real_c(1) - sleepBias )*( v_.sqrLength() + w_.sqrLength() );

   // Activating the sleep mode of the rigid body
   if( motion_ < sleepThreshold ) {
      awake_ = false;
      v_.reset();
      w_.reset();
   }
}
//*************************************************************************************************

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the rigid body currently has a supervising rigid body manager.
 *
 * \return \a true if the rigid body has a supervising manager, \a false if not.
 */
inline bool RigidBody::hasManager() const
{
   return manager_ != 0;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the supervising rigid body manager of the rigid body.
 *
 * \return The supervising rigid body manager.
 *
 * If the body is currently not supervised by an rigid body manager, the returned ManagerID is 0.
 */
inline ManagerID RigidBody::getManager()
{
   return manager_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the supervising rigid body manager of the rigid body.
 *
 * \return The supervising rigid body manager.
 *
 * If the body is currently not supervised by an rigid body manager, the returned ManagerID is 0.
 */
inline ConstManagerID RigidBody::getManager() const
{
   return manager_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body is contained in a superordinate body.
 *
 * \return \a true if the rigid body is contained in a superordinate body, \a false if not.
 */
inline bool RigidBody::hasSuperBody() const
{
   return sb_ != this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the superordinate body in which the rigid body is contained.
 *
 * \return The superordinate body.
 *
 * This function returns the direct superordinate body in which the rigid body is contained.
 * If the rigid body is not contained in another body, the returned BodyID is the body itself.
 */
inline BodyID RigidBody::getSuperBody()
{
   return sb_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the superordinate body in which the rigid body is contained.
 *
 * \return The superordinate body.
 *
 * This function returns the direct superordinate body in which the rigid body is contained.
 * If the rigid body is not contained in another body, the returned BodyID is the body itself.
 */
inline ConstBodyID RigidBody::getSuperBody() const
{
   return sb_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the top level superordinate body in which the rigid body is contained.
 *
 * \return The top level superordinate body.
 *
 * This function returns the top level superordinate body in which the rigid body is contained.
 * If the rigid body is not contained in another body, the returned BodyID is the body itself.
 */
inline BodyID RigidBody::getTopSuperBody()
{
   if( !hasSuperBody() )
      return sb_;
   else return sb_->getTopSuperBody();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the top level superordinate body in which the rigid body is contained.
 *
 * \return The superordinate body.
 *
 * This function returns the top level superordinate body in which the rigid body is contained.
 * If the rigid body is not contained in another body, the returned BodyID is the body itself.
 */
inline ConstBodyID RigidBody::getTopSuperBody() const
{
   if( !hasSuperBody() )
      return sb_;
   else return sb_->getTopSuperBody();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body is finite or not.
 *
 * \return \a true in case of a finite rigid body, \a false in case of a infinite rigid body.
 */
inline bool RigidBody::isFinite() const
{
   return finite_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body is awake or not.
 *
 * \return \a true in case the rigid body is awake, \a false if is not.
 */
inline bool RigidBody::isAwake() const
{
   if( !hasSuperBody() )
      return awake_;
   else return sb_->isAwake();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body's position is fixed or not.
 *
 * \return \a true in case of a fixed/stationary rigid body, \a false in case of a free body.
 */
inline bool RigidBody::isFixed() const
{
   return !isCommunicating() && hasInfiniteMass() && !isGlobal();
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Checks if a body is "mobile" e.g. will be integrated by the simulation
 */
inline bool RigidBody::hasInfiniteMass() const
{
   return std::isinf(getMass());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body is visible or not.
 *
 * \return \a true in case of a visible rigid body, \a false in case of an invisible body.
 */
inline bool RigidBody::isVisible() const
{
   return visible_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the unique system-specific ID of the rigid body.
 *
 * \return The system-specific ID.
 */
inline id_t RigidBody::getSystemID() const
{
   return sid_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the user-specific ID of the rigid body.
 *
 * \return The user-specific ID.
 */
inline id_t RigidBody::getID() const
{
   return uid_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the relative position of the rigid body within the superordinate body.
 *
 * \return The relative position of the rigid body.
 *
 * If the rigid body is not contained in a superordinate body, the returned relative position will
 * be \f$ \left(\begin{array}{*{3}{c}} 0 & 0 & 0 \end{array}\right) \f$.
 */
inline const Vec3& RigidBody::getRelPosition() const
{
   return rpos_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the global position of the center of mass of the rigid body.
 *
 * \return The global position of the center of mass.
 */
inline const Vec3& RigidBody::getPosition() const
{
   return gpos_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the relative linear velocity of the rigid body.
 *
 * \return The relative linear velocity.
 *
 * This function returns the linear velocity of the center of mass of the rigid body in reference
 * to the body's own frame of reference.
 */
inline const Vec3 RigidBody::getRelLinearVel() const
{
   if( hasSuperBody() )
      v_ = sb_->velFromWF( gpos_ );
   return R_.getTranspose() * v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the global linear velocity of the rigid body.
 *
 * \return The global linear velocity.
 *
 * This function returns the linear velocity of the center of mass of the rigid body in reference
 * to the global world frame.
 */
inline const Vec3& RigidBody::getLinearVel() const
{
   if( hasSuperBody() )
      v_ = sb_->velFromWF( gpos_ );
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the relative angular velocity.
 *
 * \return The relative angular velocity.
 *
 * This function returns the angluar velocity of the center of mass in reference to the body's
 * own frame of reference.
 */
inline const Vec3 RigidBody::getRelAngularVel() const
{
   if( hasSuperBody() )
      w_ = sb_->getAngularVel();
   return R_.getTranspose() * w_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the global angular velocity.
 *
 * \return The global angular velocity.
 *
 * This function returns the angluar velocity of the center of mass in reference to the global
 * world frame.
 */
inline const Vec3& RigidBody::getAngularVel() const
{
   if( hasSuperBody() )
      w_ = sb_->getAngularVel();
   return w_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the orientation of the rigid body.
 *
 * \return The orientation of the rigid body.
 *
 * This function returns the quaternion of the rigid body, which represents the orientation of
 * the body in reference to the global world frame.
 */
inline const Quat& RigidBody::getQuaternion() const
{
   return q_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the rotation of the rigid body.
 *
 * \return The rotation of the rigid body.
 *
 * This function returns the rotation matrix of the rigid body, which represents the rotation of
 * the body in reference to the global world frame.
 */
inline const Mat3& RigidBody::getRotation() const
{
   return R_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total mass of the rigid body.
 *
 * \return The total mass of the rigid body.
 */
inline real_t RigidBody::getMass() const
{
   return mass_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the inverse total mass of the rigid body.
 *
 * \return The inverse total mass of the rigid body.
 */
inline real_t RigidBody::getInvMass() const
{
   return invMass_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the volume of the rigid body.
 *
 * \return The volume of the rigid body.
 */
inline real_t RigidBody::getVolume() const
{
   WALBERLA_ABORT("getVolume is not implemented for this rigid body!");
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the moment of inertia in reference to the body frame of the rigid body.
 *
 * \return The body relative moment of inertia.
 */
inline const Mat3& RigidBody::getBodyInertia() const
{
   return I_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the moment of inertia in reference to the global world frame.
 *
 * \return The global moment of inertia.
 */
inline const Mat3 RigidBody::getInertia() const
{
   return math::transformMatrixRART(R_, I_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the inverse moment of inertia in reference to the body frame of the rigid body.
 *
 * \return The inverse body relative moment of inertia.
 */
inline const Mat3& RigidBody::getInvBodyInertia() const
{
   return Iinv_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the inverse moment of inertia in reference to the global world frame.
 *
 * \return The inverse global moment of inertia.
 */
inline const Mat3 RigidBody::getInvInertia() const
{
   return math::transformMatrixRART(R_, Iinv_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the axis-aligned bounding box of the rigid body.
 *
 * \return The axis-aligned bounding box of the rigid body.
 */
inline const AABB& RigidBody::getAABB() const
{
   return aabb_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the length of the longest side of the AABB of the rigid body.
 *
 * \return The length of the longest side of the AABB of the rigid body.
 */
inline real_t RigidBody::getAABBSize() const
{
   real_t size  = aabb_.xSize();

   if( aabb_.ySize() > size ) size = aabb_.ySize();
   if( aabb_.zSize() > size ) size = aabb_.zSize();

   return size;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the kinetic energy of the rigid body
 *
 * \return kinetic energy
 */
inline real_t   RigidBody::getKineticEnergy()      const
{
   return real_c(0.5) * getMass() * getLinearVel() * getLinearVel();
}
//*************************************************************************************************
//*************************************************************************************************
/*!\brief Returns the rotational energy of the rigid body
 *
 * \return rotational energy
 */
inline real_t   RigidBody::getRotationalEnergy()      const
{
   return real_c(0.5) * getAngularVel() * (getInertia() * getAngularVel());
}
//*************************************************************************************************
//*************************************************************************************************
/*!\brief Returns the energy of the rigid body
 *
 * \return energy. Takes into account rotational energy, kinetic energy and potential energy.
 */
inline real_t   RigidBody::getEnergy()      const
{
   return getKineticEnergy() + getRotationalEnergy();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a relative to a global vector.
 *
 * \param vx The x-component of the relative vector.
 * \param vy The y-component of the relative vector.
 * \param vz The z-component of the relative vector.
 * \return The transformed global vector.
 *
 * The function calculates the transformation of a vector in body frame to a vector in global
 * world frame.
 */
inline const Vec3 RigidBody::vectorFromBFtoWF( real_t vx, real_t vy, real_t vz ) const
{
   return R_ * Vec3( vx, vy, vz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a relative to a global vector.
 *
 * \param v the relative vector.
 * \return The transformed global vector.
 *
 * The function calculates the transformation of a vector in body frame to a vector in global
 * world frame.
 */
inline const Vec3 RigidBody::vectorFromBFtoWF( const Vec3& v ) const
{
   return R_ * v;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a relative to a global coordinate.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return The transformed global coordinate.
 *
 * The function calculates the transformation of a point relative to the body's center of
 * mass to a point in global coordinates.
 */
inline const Vec3 RigidBody::pointFromBFtoWF( real_t px, real_t py, real_t pz ) const
{
   return gpos_ + ( R_ * Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a relative to a global coordinate.
 *
 * \param rpos the relative coordinate.
 * \return The transformed global coordinate.
 *
 * The function calculates the transformation of a point relative to the body's center of
 * mass to a point in global coordinates.
 */
inline const Vec3 RigidBody::pointFromBFtoWF( const Vec3& rpos ) const
{
   return gpos_ + ( R_ * rpos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a global to a relative vector.
 *
 * \param vx The x-component of the global vector.
 * \param vy The y-component of the global vector.
 * \param vz The z-component of the global vector.
 * \return The transformed relative vector.
 *
 * The function calculates the transformation of a vector in global world frame to a vector
 * in body frame.
 */
inline const Vec3 RigidBody::vectorFromWFtoBF( real_t vx, real_t vy, real_t vz ) const
{
   return R_.getTranspose() * Vec3( vx, vy, vz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a global to a relative vector.
 *
 * \param v The global vector.
 * \return The transformed relative vector.
 *
 * The function calculates the transformation of a vector in global world frame to a vector
 * in body frame.
 */
inline const Vec3 RigidBody::vectorFromWFtoBF( const Vec3& v ) const
{
   return R_.getTranspose() * v;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a global to a relative coordinate.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return The transformed relative coordinate.
 *
 * The function calculates the transformation of a point in global coordinates to a point
 * relative to the body's center of mass.
 */
inline const Vec3 RigidBody::pointFromWFtoBF( real_t px, real_t py, real_t pz ) const
{
   return R_.getTranspose() * ( Vec3( px, py, pz ) - gpos_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transformation from a global to a relative coordinate.
 *
 * \param gpos The global coordinate.
 * \return The transformed relative coordinate.
 *
 * The function calculates the transformation of a point in global coordinates to a point
 * relative to the body's center of mass.
 */
inline const Vec3 RigidBody::pointFromWFtoBF( const Vec3& gpos ) const
{
   return R_.getTranspose() * ( gpos - gpos_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global acceleration of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return The global acceleration.
 *
 * The function calculates the global acceleration of a point in global coordinates.
 */
const Vec3 RigidBody::accFromWF( real_t px, real_t py, real_t pz ) const
{
   return accFromWF( Vec3( px, py, pz ) );
}
//*************************************************************************************************

inline id_t RigidBody::getTypeID() const
{
   WALBERLA_ASSERT_LESS( typeID_, std::numeric_limits<id_t>::max(), "You are requesting the type " \
                         " id of a body, but the static type id for the body has not yet been " \
                         " initialized! Call SetBodyTypeIDs<BodyTypeTuple>::execute to initialize!" );
   return typeID_;
}




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting the relative linear velocity of the rigid body.
 *
 * \param vx The x-component of the relative linear velocity.
 * \param vy The y-component of the relative linear velocity.
 * \param vz The z-component of the relative linear velocity.
 * \return void
 *
 * This function sets the linear velocity of the rigid body in reference to the body's own
 * frame of reference. The given relative velocity is translated into the global world frame
 * depending on the orientation of the rigid body. If the body is contained in a superordinate body
 * the function has no effect.
 *
 * In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 * on one process may invalidate the settings of the rigid body on another process. In order to
 * synchronize all rigid bodies after local changes, the simulation has to be synchronized
 * by the user. Note that any changes on remote rigid
 * bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::setRelLinearVel( real_t vx, real_t vy, real_t vz )
{
   if( !hasSuperBody() )
   {
      v_ = R_ * Vec3( vx, vy, vz );
      wake();
   }
}
//*************************************************************************************************

/// \see setRelLinearVel( real_t vx, real_t vy, real_t vz )
inline void RigidBody::setRelLinearVel( const Vec3& lvel )
{
   setRelLinearVel(lvel[0], lvel[1], lvel[2]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the global linear velocity of the rigid body.
 *
 * \param vx The x-component of the global linear velocity.
 * \param vy The y-component of the global linear velocity.
 * \param vz The z-component of the global linear velocity.
 * \return void
 *
 * This function sets the linear velocity of the rigid body. If the body is contained in a
 * superordinate body the function has no effect.
 *
 * In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 * on one process may invalidate the settings of the rigid body on another process. In order to
 * synchronize all rigid bodies after local changes, the simulation has to be synchronized
 * by the user. Note that any changes on remote rigid
 * bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::setLinearVel( real_t vx, real_t vy, real_t vz )
{
   if( !hasSuperBody() )
   {
      v_ = Vec3( vx, vy, vz );
      wake();
   }
}
//*************************************************************************************************


/// /see setLinearVel( real_t vx, real_t vy, real_t vz )
inline void RigidBody::setLinearVel( const Vec3& lvel )
{
   setLinearVel(lvel[0], lvel[1], lvel[2]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the relative angular velocity of the rigid body.
 *
 * \param ax The x-component of the relative angular velocity.
 * \param ay The y-component of the relative angular velocity.
 * \param az The z-component of the relative angular velocity.
 * \return void
 *
 * This function sets the angular velocity of the rigid body in reference to the body's own
 * frame of reference. The given relative velocity is translated into the global world frame
 * depending on the orientation of the rigid body. If the body is contained in a superordinate body
 * the function has no effect.
 *
 * In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 * on one process may invalidate the settings of the rigid body on another process. In order to
 * synchronize all rigid bodies after local changes, the simulation has to be synchronized
 * by the user. Note that any changes on remote rigid
 * bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::setRelAngularVel( real_t ax, real_t ay, real_t az )
{
   if( !hasSuperBody() )
   {
      w_ = R_ * Vec3( ax, ay, az );
      wake();
   }
}
//*************************************************************************************************


/// /see setRelAngularVel( real_t ax, real_t ay, real_t az )
inline void RigidBody::setRelAngularVel( const Vec3& avel )
{
   setRelAngularVel(avel[0], avel[1], avel[2]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the global angular velocity of the rigid body.
 *
 * \param ax The x-component of the global angular velocity.
 * \param ay The y-component of the global angular velocity.
 * \param az The z-component of the global angular velocity.
 * \return void
 *
 * This function sets the angular velocity of the rigid body. If the body is contained in a
 * superordinate body the function has no effect.
 *
 * In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 * on one process may invalidate the settings of the rigid body on another process. In order to
 * synchronize all rigid bodies after local changes, the simulation has to be synchronized
 * by the user. Note that any changes on remote rigid
 * bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::setAngularVel( real_t ax, real_t ay, real_t az )
{
   if( !hasSuperBody() )
   {
      w_ = Vec3( ax, ay, az );
      wake();
   }
}
//*************************************************************************************************


/// /see setAngularVel( real_t ax, real_t ay, real_t az )
inline void RigidBody::setAngularVel( const Vec3& avel )
{
   setAngularVel(avel[0], avel[1], avel[2]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the relative position within a superordinate body.
 *
 * \return void
 *
 * The function calculates the relative position depending on its current global position, the
 * current global position of the superordinate body and the rotation of the superordinate body.
 */
inline void RigidBody::calcRelPosition()
{
   rpos_ = sb_->R_.getTranspose()*( gpos_ - sb_->gpos_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  FORCE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the rigid body has non-zero acting forces or torques.
 *
 * \return \a true if the currently acting force and/or torque is non-zero, \a false if not.
 */
inline bool RigidBody::hasForce() const
{
   return ( force_ != real_c(0) || torque_ != real_c(0) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current force acting on the body's center of mass.
 *
 * \return The acting force.
 *
 * This function returns the current force acting on the rigid body. If the body is contained
 * in a union, then this force is part of the total force acting on the union.\n
 */
inline const Vec3& RigidBody::getForce() const
{
   return force_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current torque acting on the body's center of mass.
 *
 * \return The acting torque.
 *
 * This function returns the current torque acting in the center of mass of the rigid body. If
 * the body is contained in a union, then this torque represents the part of the total torque
 * acting on the union that results from the forces on this body.\n
 */
inline const Vec3& RigidBody::getTorque() const
{
   return torque_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Sets the super body for the current body.
 *
 *  \attention Behaviour changed compared to setSuperBody of old PE!!!
 */
inline  void RigidBody::setSB  ( BodyID body )
{
   sb_ = body;
}

//*************************************************************************************************
/*!\brief Resets the super body for the current body.
 *
 *  \attention Behaviour changed compared to setSuperBody of old PE!!!
 */
inline  void RigidBody::resetSB()
{
   sb_ = this;
}


//*************************************************************************************************
/*!\brief Set the force acting at the body's center of mass.
 *
 * \param f The acting force.
 * \return void
 */
inline void RigidBody::setForce( const Vec3& f )
{
   // Increasing the force on this rigid body
   force_ = f;

   wake();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Set the torque acting at the body's center of mass.
 *
 * \param tau The acting torque.
 * \return void
 */
inline void RigidBody::setTorque( const Vec3& tau )
{
   // Increasing the force on this rigid body
   torque_ = tau;

   wake();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the total force acting in the body's center of mass.
 *
 * \param fx The x-component of the relative force.
 * \param fy The y-component of the relative force.
 * \param fz The z-component of the relative force.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting
 * directly in the body's center of mass and increases the total force acting on the body.
 * If the body is part of a superordinate body, the force is also acting on the superordinate
 * body. Depending on the position of the superordinate body's center of mass, the force can
 * also cause a torque in the superordinate body.
 */
inline void RigidBody::addRelForce( real_t fx, real_t fy, real_t fz )
{
   addForce( R_ * Vec3( fx, fy, fz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting relative force.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting
 * directly in the body's center of mass and increases the total force acting on the body.
 * If the body is part of a superordinate body, the force is also acting on the superordinate
 * body. Depending on the position of the superordinate body's center of mass, the force can
 * also cause a torque in the superordinate body.
 */
inline void RigidBody::addRelForce( const Vec3& f )
{
   addForce( R_ * f );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param fx The x-component of the force.
 * \param fy The y-component of the force.
 * \param fz The z-component of the force.
 * \return void
 *
 * The function applies a global force to the rigid body. The given force is acting directly
 * in the body's center of mass and increases the total acting force on the body. If the rigid
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 * Depending on the position of the superordinate body's center of mass, the force can also
 * cause a torque in the superordinate body.
 */
inline void RigidBody::addForce( real_t fx, real_t fy, real_t fz )
{
   addForce( Vec3( fx, fy, fz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting force.
 * \return void
 *
 * The function applies a global force to the rigid body. The given force is acting directly
 * in the body's center of mass and increases the total acting force on the body. If the rigid
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 * Depending on the position of the superordinate body's center of mass, the force can also
 * cause a torque in the superordinate body.
 */
inline void RigidBody::addForce( const Vec3& f )
{
   if (sb_ == this)
   {
      // Increasing the force on this rigid body
      force_ += f;
      wake();
   } else
   {
      sb_->addForceAtPos( f, gpos_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param fx The x-component of the relative force.
 * \param fy The y-component of the relative force.
 * \param fz The z-component of the relative force.
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting at
 * the specified body relative coordinate and increases the total acting force on the rigid
 * body. Depending on the position, the force can cause a torque in the body's center of mass.
 * If the body is part of a superordinate body, the force is also acting on the superordinate
 * body.
 */
inline void RigidBody::addRelForceAtRelPos( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz )
{
   addForceAtPos( vectorFromBFtoWF( fx, fy, fz ), pointFromBFtoWF( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting relative force.
 * \param rpos The relative coordinate.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting at
 * the specified body relative coordinate and increases the total acting force on the rigid
 * body. Depending on the position, the force can cause a torque in the body's center of mass.
 * If the body is part of a superordinate body, the force is also acting on the superordinate
 * body.
 */
inline void RigidBody::addRelForceAtRelPos( const Vec3& f, const Vec3& rpos )
{
   addForceAtPos( vectorFromBFtoWF( f ), pointFromBFtoWF( rpos ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param fx The x-component of the relative force.
 * \param fy The y-component of the relative force.
 * \param fz The z-component of the relative force.
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting at
 * the specified global coordinate and increases the total acting force on the rigid body.
 * Depending on the position, the force can cause a torque in the body's center of mass. If the
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 */
inline void RigidBody::addRelForceAtPos( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz )
{
   addForceAtPos( vectorFromBFtoWF( fx, fy, fz ), Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting relative force.
 * \param gpos The global coordinate.
 * \return void
 *
 * The function applies a body relative force to the rigid body. The given force is acting at
 * the specified global coordinate and increases the total acting force on the rigid body.
 * Depending on the position, the force can cause a torque in the body's center of mass. If the
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 */
inline void RigidBody::addRelForceAtPos( const Vec3& f, const Vec3& gpos )
{
   addForceAtPos( vectorFromBFtoWF( f ), gpos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param fx The x-component of the global force.
 * \param fy The y-component of the global force.
 * \param fz The z-component of the global force.
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return void
 *
 * The function applies a global force to the rigid body. The given force is acting at the
 * specified body relative coordinate and increases the total acting force on the rigid body.
 * Depending on the position, the force can cause a torque in the body's center of mass. If the
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 */
inline void RigidBody::addForceAtRelPos( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz )
{
   addForceAtPos( Vec3( fx, fy, fz ), pointFromBFtoWF( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting global force.
 * \param rpos The relative coordinate.
 * \return void
 *
 * The function applies a global force to the rigid body. The given force is acting at the
 * specified body relative coordinate and increases the total acting force on the rigid body.
 * Depending on the position, the force can cause a torque in the body's center of mass. If the
 * body is part of a superordinate body, the force is also acting on the superordinate body.
 */
inline void RigidBody::addForceAtRelPos( const Vec3& f, const Vec3& rpos )
{
   addForceAtPos( f, pointFromBFtoWF( rpos ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param fx The x-component of the global force.
 * \param fy The y-component of the global force.
 * \param fz The z-component of the global force.
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return void
 *
 * The given force is acting at the specified coordinate and increases the total acting
 * force on the rigid body. Depending on the position, the force can cause a torque in the
 * body's center of mass. If the body is part of a superordinate body, the force is also
 * acting on the superordinate body.
 */
inline void RigidBody::addForceAtPos( real_t fx, real_t fy, real_t fz, real_t px, real_t py, real_t pz )
{
   addForceAtPos( Vec3( fx, fy, fz ), Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the force acting in the body's center of mass.
 *
 * \param f The acting global force.
 * \param gpos The global coordinate.
 * \return void
 *
 * The given force is acting at the specified coordinate and increases the total acting
 * force on the rigid body. Depending on the position, the force can cause a torque in the
 * body's center of mass. If the body is part of a superordinate body, the force is also
 * acting on the superordinate body.
 */
inline void RigidBody::addForceAtPos( const Vec3& f, const Vec3& gpos )
{
   if (sb_ == this)
   {
      // Increasing the force and torque on this rigid body
      force_  += f;
      torque_ += ( gpos - gpos_ ) % f;
      wake();
   } else
   {
      sb_->addForceAtPos( f, gpos );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increasing the torque acting in the body's center of mass.
 *
 * \param tx The x-component of the torque.
 * \param ty The y-component of the torque.
 * \param tz The z-component of the torque.
 * \return void
 *
 * The torque is acting directly in the body's center of mass and increases the total acting
 * torque on the body. If the rigid body is part of a superordinate body, the torque is applied
 * to the superordinate body instead. It is \b not possible to apply a torque on subordinate
 * rigid bodies!
 */
inline void RigidBody::addTorque( real_t tx, real_t ty, real_t tz )
{
   addTorque( Vec3( tx, ty, tz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increasing the torque acting in the body's center of mass.
 *
 * \param t The acting torque.
 * \return void
 *
 * The torque is acting directly in the body's center of mass and increases the total acting
 * torque on the body. If the rigid body is part of a superordinate body, the torque is applied
 * to the superordinate body instead. It is \b not possible to apply a torque on subordinate
 * rigid bodies!
 */
inline void RigidBody::addTorque( const Vec3& t )
{
   if( !hasSuperBody() ) {
      torque_ += t;
      wake();
   }
   else sb_->addTorque( t );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting all acting forces and torques from the rigid body.
 *
 * \return void
 */
inline void RigidBody::resetForceAndTorque()
{
   force_  = Vec3(0);
   torque_ = Vec3(0);
}
//*************************************************************************************************




//=================================================================================================
//
//  IMPULSE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Applying an impulse in the body's center of mass.
 *
 * \param jx The x-component of the impulse.
 * \param jy The y-component of the impulse.
 * \param jz The z-component of the impulse.
 * \return void
 *
 * The impulse is acting directly in the body's center of mass and instantaneously changes
 * the linear velocity of the rigid body. If the body is part of a superordinate body, the
 * impulse is also acting on the superordinate body. Depending on the position of the
 * superordinate body's center of mass, the impulse can also change the angular velocity
 * of the rigid body (and the superordinate body).
 */
inline void RigidBody::addImpulse( real_t jx, real_t jy, real_t jz )
{
   addImpulse( Vec3( jx, jy, jz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applying an impulse in the body's center of mass.
 *
 * \param j The acting impulse.
 * \return void
 *
 * The impulse is acting directly in the body's center of mass and instantaneously changes
 * the linear velocity of the rigid body. If the body is part of a superordinate body, the
 * impulse is also acting on the superordinate body. Depending on the position of the
 * superordinate body's center of mass, the impulse can also change the angular velocity
 * of the rigid body (and the superordinate body).
 */
inline void RigidBody::addImpulse( const Vec3& j )
{
   if( !hasSuperBody() ) {
      v_ += j * invMass_;
      wake();
   }
   else sb_->addImpulseAtPos( j, gpos_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applying an impulse at the given global coordinate.
 *
 * \param jx The x-component of the impulse.
 * \param jy The y-component of the impulse.
 * \param jz The z-component of the impulse.
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return void
 *
 * The given impulse is acting at the specified coordinate and instantaneously changes the linear
 * velocity of the rigid body. Depending on the position of the body's center of mass, the impulse
 * can also change the angular velocity. If the rigid body is part of a superordinate body, the
 * impulse is also acting on the superordinate body.
 */
inline void RigidBody::addImpulseAtPos( real_t jx, real_t jy, real_t jz, real_t px, real_t py, real_t pz )
{
   addImpulseAtPos( Vec3( jx, jy, jz ), Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applying an impulse at the given global coordinate.
 *
 * \param j The acting impulse.
 * \param p The global coordinate.
 * \return void
 *
 * The given impulse is acting at the specified coordinate and instantaneously changes the linear
 * velocity of the rigid body. Depending on the position of the body's center of mass, the impulse
 * can also change the angular velocity. If the rigid body is part of a superordinate body, the
 * impulse is also acting on the superordinate body.
 */
inline void RigidBody::addImpulseAtPos( const Vec3& j, const Vec3& p )
{
   if( !hasSuperBody() ) {
      v_ += j * invMass_;
      w_ += getInvInertia() * ( ( p - gpos_ ) % j );
      wake();
   }
   else sb_->addImpulseAtPos( j, p );
}
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR INTERNAL CHANGES IN COMPOUND GEOMETRIES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Signals an internal modification of a contained subordiante body.
 *
 * \return void
 *
 * This function can be used by derived primitive geometries to signal an internal modification
 * to its superordinate body.
 */
inline void RigidBody::signalModification()
{
   if( hasSuperBody() ) sb_->handleModification();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Signals a position change of a contained subordiante body.
 *
 * \return void
 *
 * This function can be used by derived primitive geometries to signal a translation or
 * a position change to its superordinate body.
 */
inline void RigidBody::signalTranslation()
{
   if( hasSuperBody() ) sb_->handleTranslation();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Signals an orientation change of a contained subordiante body.
 *
 * \return void
 *
 * This function can be used by derived primitive geometries to signal a rotation or
 * orientation change to its superordinate body.
 */
inline void RigidBody::signalRotation()
{
   if( hasSuperBody() ) sb_->handleRotation();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Signals a fixation change of a contained subordiante body.
 *
 * \return void
 *
 * This function can be used by derived primitive geometries to signal a change of the fixation
 * flag to its superordinate body.
 */
inline void RigidBody::signalFixation()
{
   //if( hasSuperBody() && sb_->isFixed() != isFixed() ) sb_->handleFixation();
}
//*************************************************************************************************




//=================================================================================================
//
//  MPI FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the rigid body is remote or not.
 *
 * \return \a true in case of a remote rigid body, \a false in case of a local body.
 */
inline bool RigidBody::isRemote() const
{
   return remote_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the remote flag of the rigid body.
 *
 * \param remote \a true to declare the rigid body remote, \a false declare it local.
 * \return void
 *
 * This function sets the remote flag of the rigid body. Note that this function should not be
 * used explicitly, but is automatically called during the MPI communication to set the remote
 * status of a rigid body within the simulation world. Using this function explicitly may lead
 * to simulation errors during a parallel simulation!
 */
inline void RigidBody::setRemote( bool remote )
{
   remote_ = remote;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the rigid body is local or not.
 *
 * \return \a true in case of a local rigid body, \a false otherwise.
 */
inline bool RigidBody::isCommunicating() const
{
   return communicating_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the local flag of the rigid body.
 *
 * \return void
 *
 * This function declares the rigid body as a local body. Note that this function should not be
 * used explicitly, but is automatically called during the setup of local rigid bodies. Using
 * this function explicitly may lead to simulation errors during a parallel simulation!
 */
inline void RigidBody::setCommunicating( const bool communicating )
{
   communicating_ = communicating;    // Marking the rigid body as local
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns whether the rigid body is global or not.
 *
 * \return \a true in case of a global rigid body, \a false otherwise.
 */
inline bool RigidBody::isGlobal() const
{
   return global_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the global flag of the rigid body.
 *
 * \return void
 *
 * This function declares the rigid body as a global body. Note that this function should not be
 * used explicitly, but is automatically called during the setup of local rigid bodies. Using
 * this function explicitly may lead to simulation errors during a parallel simulation!
 */
inline void RigidBody::setGlobal(const bool global)
{
   global_ = global;    // Marking the rigid body as global
}
//*************************************************************************************************



//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the global acceleration of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return The global acceleration.
 *
 * The function calculates the global acceleration of a point in global coordinates.
 */
inline const Vec3 RigidBody::accFromWF( const Vec3& gpos ) const
{
   if( !hasSuperBody() )
   {
      // Calculating the linear acceleration
      const Vec3 vdot( force_ * invMass_ );

      // Calculating the angular acceleration
      const Vec3 wdot( getInvInertia() * ( torque_ - w_ % ( getInertia() * w_ ) ) );

      // Calculating the distance to the center of mass of the superordinate body
      const Vec3 r( gpos - gpos_ );

      // Calculating the acceleration of the point 'gpos'
      return vdot + wdot % r + w_ % ( w_ % r );
   }
   else return sb_->accFromWF( gpos );
}
//*************************************************************************************************




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================

inline void RigidBody::setFinite( const bool finite )
{
   finite_ = finite;
}


//*************************************************************************************************
/*!\brief Sets mass and inertia of a rigid body. Also calculates inverse values.
 *
 * \param mass mass to be set (may be infinity)
 * \param inertia inertia to be set (if mass is infinity this one will be ignored)
 */
inline void RigidBody::setMassAndInertia( const real_t mass, const Mat3& inertia )
{
   if ( std::isinf( mass ) )
   {
      // Adjusting the inverse mass and inverse moment of inertia
      mass_    = std::numeric_limits<real_t>::infinity();
      I_       = Mat3::makeDiagonalMatrix(std::numeric_limits<real_t>::infinity());
      invMass_ = real_c(0);
      Iinv_    = Mat3( real_c(0) );
   } else
   {
      // Adjusting the inverse mass and inverse moment of inertia
      mass_    = mass;
      I_       = inertia;
      invMass_ = real_c(1) / mass_;
      Iinv_    = I_.getInverse();
   }

   if (hasSuperBody()) WALBERLA_LOG_WARNING("Changing the mass of a body contained in a Union is currently untested!!!");
   // Signaling the fixation change to the superordinate body
   signalFixation();
}
//*************************************************************************************************



//*************************************************************************************************
/*!\fn void RigidBody::setVisible( bool visible )
 * \brief Setting the rigid body visible/invisible.
 *
 * \param visible \a true to make the rigid body visible, \a false to make it invisible.
 * \return void
 */
//*************************************************************************************************
inline void RigidBody::setVisible( bool visible )
{
   visible_ = visible;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\fn void RigidBody::setPosition( real_t px, real_t py, real_t pz )
 * \brief Setting the global position of the rigid body.
 *
 * \param px The x-component of the global position.
 * \param py The y-component of the global position.
 * \param pz The z-component of the global position.
 * \return void
 *
 * This function sets the global position of the rigid body to the given coordinate (px,py,pz).
 *
 * \b Note:
 * - Setting the position of a rigid body contained in a union changes the mass distribution and
 *   geometry of the union. Therefore this may cause an invalidation of links contained in the
 *   union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the body on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setPositionImpl( real_t px, real_t py, real_t pz )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to move a fixed body: " << *this);

   gpos_[0] = px;
   gpos_[1] = py;
   gpos_[2] = pz;

   calcBoundingBox();    // Updating the axis-aligned bounding box of the box
   wake();               // Waking the box from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\fn void RigidBody::setPosition( real_t px, real_t py, real_t pz )
 * \brief Setting the global position of the rigid body.
 *
 * \param px The x-component of the global position.
 * \param py The y-component of the global position.
 * \param pz The z-component of the global position.
 * \return void
 *
 * This function sets the global position of the rigid body to the given coordinate (px,py,pz).
 *
 * \b Note:
 * - Setting the position of a rigid body contained in a union changes the mass distribution and
 *   geometry of the union. Therefore this may cause an invalidation of links contained in the
 *   union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the body on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setPosition( real_t px, real_t py, real_t pz )
{
   setPositionImpl(px, py, pz);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\fn void RigidBody::setPosition( const Vec3& gpos )
 * \brief Setting the global position of the rigid body.
 *
 * \param gpos The global position.
 * \return void
 *
 * This function sets the global position of the rigid body to the given coordinate \a gpos.
 *
 * \b Note:
 * - Setting the position of a rigid body contained in a union changes the mass distribution and
 *   geometry of the union. Therefore this may cause an invalidation of links contained in the
 *   union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the body on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setPosition( const Vec3& gpos )
{
   setPositionImpl( gpos[0], gpos[1], gpos[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\fn void RigidBody::setOrientationImpl( real_t r, real_t i, real_t j, real_t k )
 * \brief Setting the global orientation of the rigid body.
 *
 * \param r The value for the real_t part.
 * \param i The value for the first imaginary part.
 * \param j The value for the second imaginary part.
 * \param k The value for the third imaginary part.
 * \return void
 *
 * This function sets the global orientation of the rigid body to the given quaternion
 * \a (r,i,j,k).
 *
 * \b Note:
 * - Setting the orientation of a rigid body contained in a union changes the mass distribution
 *   and geometry of the union. Therefore this changes the union and may cause an invalidation
 *   of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the box on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setOrientationImpl( real_t r, real_t i, real_t j, real_t k )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to rotate a fixed body: " << *this);

   q_.set(r, i, j, k);          // Setting the orientation of the box
   R_ = q_.toRotationMatrix();  // Updating the rotation of the box

   calcBoundingBox();  // Updating the axis-aligned bounding box of the box
   wake();             // Waking the box from sleep mode
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\fn void RigidBody::setOrientation( real_t r, real_t i, real_t j, real_t k )
 * \brief Setting the global orientation of the rigid body.
 *
 * \param r The value for the real_t part.
 * \param i The value for the first imaginary part.
 * \param j The value for the second imaginary part.
 * \param k The value for the third imaginary part.
 * \return void
 *
 * This function sets the global orientation of the rigid body to the given quaternion
 * \a (r,i,j,k).
 *
 * \b Note:
 * - Setting the orientation of a rigid body contained in a union changes the mass distribution
 *   and geometry of the union. Therefore this changes the union and may cause an invalidation
 *   of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the box on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setOrientation( real_t r, real_t i, real_t j, real_t k )
{
   setOrientationImpl(r, i, j, k);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\fn void RigidBody::setOrientation( const Quat& q )
 * \brief Setting the global orientation of the rigid body.
 *
 * \param q The global orientation.
 * \return void
 *
 * This function sets the global orientation of the rigid body to the given quaternion \a q.
 *
 * \b Note:
 * - Setting the orientation of a rigid body contained in a union changes the mass distribution
 *   and geometry of the union. Therefore this changes the union and may cause an invalidation
 *   of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 *   on one process may invalidate the settings of the box on another process. In order to
 *   synchronize all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid
 *   bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
//*************************************************************************************************
inline void RigidBody::setOrientation( const Quat& q )
{
   setOrientationImpl( q[0], q[1], q[2], q[3] );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\fn void RigidBody::setMassAndInertiaToInfinity( )
 * \brief Setting the mass to infinity. This will also make the inertia tensor infinite.
 *
 * \attention This cannot be undone!
 */
//*************************************************************************************************
inline void RigidBody::setMassAndInertiaToInfinity()
{
   setMassAndInertia( std::numeric_limits<real_t>::infinity(), Mat3(std::numeric_limits<real_t>::infinity()) );
}
//*************************************************************************************************



//=================================================================================================
//
//  TRANSLATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Translation of the center of mass of the rigid body by the displacement vector
 * \brief (\a dx,\a dy,\a dz).
 *
 * \param dx The x-component of the translation/displacement.
 * \param dy The y-component of the translation/displacement.
 * \param dz The z-component of the translation/displacement.
 * \return void
 *
 * \b Note:
 * - Translating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this may cause an invalidation of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::translateImpl( real_t dx, real_t dy, real_t dz )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to translate a fixed body: " << *this);

   gpos_[0] += dx;
   gpos_[1] += dy;
   gpos_[2] += dz;

   calcBoundingBox();    // Updating the axis-aligned bounding box
   wake();               // Waking the rigid body from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Translation of the center of mass of the rigid body by the displacement vector
 * \brief (\a dx,\a dy,\a dz).
 *
 * \param dx The x-component of the translation/displacement.
 * \param dy The y-component of the translation/displacement.
 * \param dz The z-component of the translation/displacement.
 * \return void
 *
 * \b Note:
 * - Translating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this may cause an invalidation of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::translate( real_t dx, real_t dy, real_t dz )
{
   translateImpl( dx, dy, dz);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Translation of the center of mass of the rigid body by the displacement vector \a dp.
 *
 * \param dp The displacement vector.
 * \return void
 *
 * \b Note:
 * - Translating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this may cause an invalidation of links contained in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::translate( const Vec3& dp )
{
   translateImpl( dp[0], dp[1], dp[2] );
}
//*************************************************************************************************


//=================================================================================================
//
//  ROTATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation of the rigid body around the global rotation axis (x,y,z) by the rotation angle \a angle.
 *
 * \param x The x-component of the global rotation axis.
 * \param y The y-component of the global rotation axis.
 * \param z The z-component of the global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * Changing the orientation/rotation of the rigid body. The rigid body is rotated around its center of mass
 * around the given axis \a (x,y,z) by \a angle degrees (radian measure).\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) box on one
 *   process may invalidate the settings of the box on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotate( real_t x, real_t y, real_t z, real_t angle )
{
   rotateImpl( Quat( Vec3(x, y, z), angle) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around the specified global rotation axis by the rotation
 * \brief angle \a angle.
 *
 * \param axis The global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * Changing the orientation/rotation of the rigid body. The rigid body is rotated around its center of mass
 * around the given axis \a axis by \a angle degrees (radian measure).\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) box on one
 *   process may invalidate the settings of the box on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotate( const Vec3& axis, real_t angle )
{
   rotateImpl( Quat( axis, angle ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body by the Euler angles \a xangle, \a yangle and \a zangle.
 *
 * \param xangle Rotation around the x-axis (radian measure).
 * \param yangle Rotation around the y-axis (radian measure).
 * \param zangle Rotation around the z-axis (radian measure).
 * \return void
 *
 * Changing the orientation/rotation of the rigid body. The rigid body is rotated around its center of mass
 * by the Euler angles \a xangle, \a yangle and \a zangle (all in radian measure). The rotations
 * are applied in the order x, y, and z.\n
 *
 * \b Note:
 * - Rotating a box contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) box on one
 *   process may invalidate the settings of the box on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotate( real_t xangle, real_t yangle, real_t zangle )
{
   Quat q;
   q.rotateX( xangle );  // Rotation around the x-axis
   q.rotateY( yangle );  // Rotation around the y-axis
   q.rotateZ( zangle );  // Rotation around the z-axis
   rotateImpl( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body by the Euler angles \a euler.
 *
 * \param euler 3-dimensional vector of the three rotation angles (radian measure).
 * \return void
 *
 * Changing the orientation/rotation of the rigid body. The rigid body is rotated around its center of mass
 * by the Euler angles \a euler (all components in radian measure). The rotations are applied
 * in the order x, y, and z.\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) box on one
 *   process may invalidate the settings of the box on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotate( const Vec3& euler )
{
   rotateImpl( Quat(euler) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body by the quaternion \a dq.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 *
 * Changing the orientation/rotation of the rigid body. The rigid body is rotated around its center of mass
 * by the quaternion \a dq. \n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) box on one
 *   process may invalidate the settings of the box on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotate( const Quat& dq )
{
   rotateImpl( dq );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Implements the rotation of a rigid body. May be overwritten in dervied classes for performance reasons.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 */
inline void RigidBody::rotateImpl( const Quat& dq )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to rotate a fixed body: " << *this);

   q_ = dq * q_;                // Updating the orientation of the box
   R_ = q_.toRotationMatrix();  // Updating the rotation of the box

   calcBoundingBox();  // Updating the axis-aligned bounding box of the box
   wake();             // Waking the box from sleep mode
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around the origin of the global world frame.
 *
 * \param x The x-component of the global rotation axis.
 * \param y The y-component of the global rotation axis.
 * \param z The z-component of the global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the origin of the global world frame and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * around the given axis \a (x,y,z) by \a angle degrees (radian measure).\n
 *
 * \see rotateAroundOrigin( const Quat& dq )
 */
inline void RigidBody::rotateAroundOrigin( real_t x, real_t y, real_t z, real_t angle )
{
   rotateAroundOrigin( Quat(Vec3( x, y, z ), angle) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around the origin of the global world frame.
 *
 * \param axis The global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the origin of the global world frame and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * around the given axis \a axis by \a angle degrees (radian measure).\n
 *
 * \see rotateAroundOrigin( const Quat& dq )
 */
inline void RigidBody::rotateAroundOrigin( const Vec3& axis, real_t angle )
{
   rotateAroundOrigin( Quat(axis, angle) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around the origin of the global world frame.
 *
 * \param xangle Rotation around the x-axis (radian measure).
 * \param yangle Rotation around the y-axis (radian measure).
 * \param zangle Rotation around the z-axis (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the origin of the global world frame and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * by the Euler angles \a xangle, \a yangle and \a zangle (all components in radian measure).
 * The rotations are applied in the order x, y, and z.\n
 *
 * \see rotateAroundOrigin( const Quat& dq )
 */
inline void RigidBody::rotateAroundOrigin( real_t xangle, real_t yangle, real_t zangle )
{
   rotateAroundOrigin( Quat(xangle, yangle, zangle) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around the origin of the global world frame.
 *
 * \param euler 3-dimensional vector of the three rotation angles (radian measure).
 * \return void
 *
 * This function rotates the box around the origin of the global world frame and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * by the Euler angles \a euler (all components in radian measure). The rotations are
 * applied in the order x, y, and z.\n
 *
 * \see rotateAroundOrigin( const Quat& dq )
 */
inline void RigidBody::rotateAroundOrigin( const Vec3& euler )
{
   rotateAroundOrigin( Quat(euler[0], euler[1], euler[2]) );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Rotation of the rigid body around the origin of the global world frame.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 *
 * This function rotates the rigid body around the origin of the global world frame and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * by the Euler angles \a euler (all components in radian measure). The rotations are
 * applied in the order x, y, and z.\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotateAroundOrigin( const Quat& dq )
{
   rotateAroundOriginImpl( dq );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Implements the rotation of a rigid body. May be overwritten in dervied classes for performance reasons.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 */
inline void RigidBody::rotateAroundOriginImpl( const Quat& dq )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to rotate a fixed body: " << *this);

   gpos_ = dq.rotate( gpos_ );     // Updating the global position of the box
   q_    = dq * q_;                // Updating the orientation of the box
   R_    = q_.toRotationMatrix();  // Updating the rotation of the box

   calcBoundingBox();    // Updating the axis-aligned bounding box of the box
   wake();               // Waking the box from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Rotation of the rigid body around a specific global coordinate.
 *
 * \param point The global center of the rotation.
 * \param axis The global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the given global coordiante \a point and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * around the given axis \a axis by \a angle degrees (radian measure).\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotateAroundPointImpl( const Vec3& point, const Quat& dq )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to rotate a fixed body: " << *this);

   const Vec3 dp( gpos_ - point );

   gpos_ = point + dq.rotate( dp );  // Updating the global position of the box
   q_    = dq * q_;                  // Updating the orientation of the box
   R_    = q_.toRotationMatrix();    // Updating the rotation of the box

   calcBoundingBox();    // Updating the axis-aligned bounding box of the box
   wake();               // Waking the box from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around a specific global coordinate.
 *
 * \param point The global center of the rotation.
 * \param axis The global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the given global coordiante \a point and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * around the given axis \a axis by \a angle degrees (radian measure).\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotateAroundPoint( const Vec3& point, const Vec3& axis, real_t angle )
{
   rotateAroundPointImpl(point, Quat(axis, angle) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the rigid body around a specific global coordinate.
 *
 * \param point The global center of the rotation.
 * \param euler 3-dimensional vector of the three rotation angles (radian measure).
 * \return void
 *
 * This function rotates the rigid body around the given global coordinate \a point and changes
 * both the global position and the orientation/rotation of the rigid body. The rigid body is rotated
 * by the Euler angles \a euler (all components in radian measure). The rotations are
 * applied in the order x, y, and z.\n
 *
 * \b Note:
 * - Rotating a rigid body contained in a union changes the mass distribution and geometry of the
 *   union. Therefore this changes the union and may cause an invalidation of links contained
 *   in the union.
 * - In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body on one
 *   process may invalidate the settings of the rigid body on another process. In order to synchronize
 *   all rigid bodies after local changes, the simulation has to be synchronized
 *   by the user. Note that any changes on remote rigid bodies are
 *   neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::rotateAroundPoint( const Vec3& point, const Vec3& euler )
{
   rotateAroundPointImpl(point, Quat(euler) );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the rigid body.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the rigid body, \a false if not.
 */
inline bool RigidBody::containsRelPoint( real_t px, real_t py, real_t pz ) const
{
   return containsRelPointImpl( px, py, pz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the rigid body.
 *
 * \param rpos The relative coordinate.
 * \return \a true if the point lies inside the rigid body, \a false if not.
 */
inline bool RigidBody::containsRelPoint( const Vec3& rpos ) const
{
   return containsRelPointImpl( rpos[0], rpos[1], rpos[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in global coordinates lies inside the rigid body.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return \a true if the point lies inside the rigid body, \a false if not.
 */
inline bool RigidBody::containsPoint( real_t px, real_t py, real_t pz ) const
{
   const Vec3& temp = pointFromWFtoBF( px, py, pz );
   return containsRelPointImpl( temp[0], temp[1], temp[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in global coordinates lies inside the rigid body.
 *
 * \param gpos The global coordinate.
 * \return \a true if the point lies inside the rigid body, \a false if not.
 */
inline bool RigidBody::containsPoint( const Vec3& gpos ) const
{
   const Vec3& temp = pointFromWFtoBF( gpos );
   return containsRelPointImpl( temp[0], temp[1], temp[2] );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the rigid body.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the rigid body, \a false if not.
 */
inline bool RigidBody::containsRelPointImpl( real_t /*px*/, real_t /*py*/, real_t /*pz*/ ) const
{
   WALBERLA_ABORT( "Contains point function not implemented!" );
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in relative coordinates lies on the surface of the rigid body.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the rigid body, \a false if not.
 *
 * The tolerance level of the check is surfaceThreshold.
 */
inline bool RigidBody::isSurfaceRelPoint( real_t px, real_t py, real_t pz ) const
{
   return isSurfaceRelPointImpl( px, py, pz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the rigid body.
 *
 * \param rpos The relative coordinate.
 * \return \a true if the point lies on the surface of the rigid body, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
inline bool RigidBody::isSurfaceRelPoint( const Vec3& rpos ) const
{
   return isSurfaceRelPointImpl( rpos[0], rpos[1], rpos[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in global coordinates lies on the surface of the rigid body.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return \a true if the point lies on the surface of the rigid body, \a false if not.
 *
 * The tolerance level of the check is surfaceThreshold.
 */
inline bool RigidBody::isSurfacePoint( real_t px, real_t py, real_t pz ) const
{
   const Vec3& temp = pointFromWFtoBF( px, py, pz );
   return isSurfaceRelPointImpl( temp[0], temp[1], temp[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in global coordinates lies on the surface of the rigid body.
 *
 * \param gpos The global coordinate.
 * \return \a true if the point lies on the surface of the rigid body, \a false if not.
 *
 * The tolerance level of the check is surfaceThreshold.
 */
inline bool RigidBody::isSurfacePoint( const Vec3& gpos ) const
{
   const Vec3& temp = pointFromWFtoBF( gpos );
   return isSurfaceRelPointImpl( temp[0], temp[1], temp[2] );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Checks, whether a point in relative coordinates lies on the surface of the rigid body.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the rigid body, \a false if not.
 *
 * The tolerance level of the check is surfaceThreshold.
 */
inline bool RigidBody::isSurfaceRelPointImpl( real_t /*px*/, real_t /*py*/, real_t /*pz*/ ) const
{
   WALBERLA_ABORT( "Surface point calculation not implemented!" );
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction a\ d.
 */
inline Vec3 RigidBody::support(const Vec3& /*d*/) const
{
   WALBERLA_ABORT( "Undefined support point calculation: \n" << this );
   return gpos_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates
 * \return The support point in world-frame coordinates in direction a\ d extended by a vector in
 *         direction \a d of length \a pe::contactThreshold.
 */
inline Vec3 RigidBody::supportContactThreshold(const Vec3& /*d*/) const
{
   WALBERLA_ABORT( "Undefined support point calculation: \n" << this  );
   return gpos_;
}
//*************************************************************************************************


//=================================================================================================
//
//  FIXATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting the global position (the center of mass) of the rigid body fixed.
 *
 * \return void
 *
 * This function fixes the global position (the center of mass) of a finite rigid body. If
 * the body is contained in a superordinate body, fixing the contained body will also fix the
 * global position of the superordinate body. In case the body is infinite or contained in an
 * infinite superordinate body (as for instance a plane or an union containing a plane) the
 * function has no effect.
 *
 * In case of a <b>MPI parallel simulation</b>, changing the settings of a (local) rigid body
 * on one process may invalidate the settings of the rigid body on another process. In order to
 * synchronize all rigid bodies after local changes, the simulation has to be synchronized
 * by the user. Note that any changes on remote rigid
 * bodies are neglected and overwritten by the settings of the rigid body on its local process!
 */
inline void RigidBody::fix()
{
   setCommunicating( false );
   setMassAndInertia( getMass(), getInertia() );

   // Signaling the fixation change to the superordinate body
   signalFixation();
}
//*************************************************************************************************


//=================================================================================================
//
//  SIMULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Translation update of a subordinate rigid body.
 *
 * \param dp Change in the global position of the superordinate rigid body.
 * \return void
 *
 * This update function is triggered by the superordinate body in case of a translational
 * movement. This movement involves a change in the global position and the axis-aligned
 * bounding box.
 */
inline void RigidBody::update( const Vec3& dp )
{
   // Checking the state of the sphere
   WALBERLA_ASSERT( checkInvariants(), "Invalid sphere state detected" );
   WALBERLA_ASSERT( hasSuperBody(), "Invalid superordinate body detected" );

   // Updating the global position
   gpos_ += dp;

   // Setting the axis-aligned bounding box
   calcBoundingBox();

   // Checking the state of the sphere
   WALBERLA_ASSERT( checkInvariants(), "Invalid sphere state detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation update of a subordinate rigid body.
 *
 * \param dq Change in the orientation of the superordinate rigid body.
 * \return void
 *
 * This update function is triggered by the superordinate body in case of a rotational movement.
 * This movement involves a change in the global position, the orientation/rotation and the
 * axis-aligned bounding box of the sphere.
 */
inline void RigidBody::update( const Quat& dq )
{
   // Checking the state of the sphere
   WALBERLA_ASSERT( checkInvariants(), "Invalid sphere state detected" );
   WALBERLA_ASSERT( hasSuperBody(), "Invalid superordinate body detected" );

   // Calculating the new global position
   gpos_ = sb_->getPosition() + ( sb_->getRotation() * rpos_ );

   // Calculating the new orientation and rotation
   q_ = dq * q_;
   R_ = q_.toRotationMatrix();

   // Setting the axis-aligned bounding box
   calcBoundingBox();

   // Checking the state of the sphere
   WALBERLA_ASSERT( checkInvariants(), "Invalid sphere state detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  SIMULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn void RigidBody::update( const Vec3& dp )
 * \brief Translation update of a subordinate rigid body.
 *
 * \param dp Change in the global position of the superordinate rigid body.
 * \return void
 *
 * This function calculates the necessary updates for a subordinate rigid body contained
 * in a superordinate rigid body that has performed a translational movement. This function
 * is triggered automatically by the superordinate body in case of a translational movement.
 * All classes deriving from the RigidBody class have to implement this function to update
 * the properties of the rigid body. The following properties of the body might change due
 * to this translation. All derived classes have to make sure these properties are updated
 * correctly:
 *
 *   - the global position
 *   - the axis-aligned bounding box
 */
//*************************************************************************************************


//*************************************************************************************************
/*!\fn void RigidBody::update( const Quat& dq )
 * \brief Rotation update of a subordinate rigid body.
 *
 * \param dq Change in the orientation of the superordinate rigid body.
 * \return void
 *
 * This function calculates the necessary updates for a subordinate rigid body contained
 * in a superordinate rigid body that has performed a rotational movement. The function
 * is triggered automatically by the superordinate body in case of a rotational movement.
 * All classes deriving from the RigidBody class have to implement this function to update
 * the properties of the rigid body. The following properties of the body might change due
 * to this rotation. All derived classes have to make sure these properties are updated
 * correctly:
 *
 *   - the global position
 *   - the orientation/rotation (i.e. the quaterion and the rotation matrix)
 *   - the axis-aligned bounding box
 */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR INTERNAL CHANGES IN COMPOUND GEOMETRIES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Handling an internal modification of a contained subordinate body.
 *
 * \return void
 *
 * This function handles an internal modification of one of the contained subordinate bodies.
 * Derived compound geometries that contain other primitive or compound bodies are required
 * to override this function in order to react to the modification. All primitive geometries
 * can use the empty default implementation.
 */
inline void RigidBody::handleModification()
{
   WALBERLA_ABORT( "Invalid call of default handle function" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Handling a position change of a contained subordinate body.
 *
 * \return void
 *
 * This function handles a translation or position change of one of the contained subordinate
 * bodies. Derived compound geometries that contain other primitive bodies are required to
 * override this function in order to react to the translation. All primitive geometries can
 * use the empty default implementation.
 */
inline void RigidBody::handleTranslation()
{
   WALBERLA_ABORT( "Invalid call of default handle function" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Handling an orientation change of a contained subordinate body.
 *
 * \return void
 *
 * This function handles a rotation or orientation change of one of the contained subordinate
 * bodies. Derived compound geometries that contain other primitive bodies are required to
 * override this function in order to react to the rotation. All primitive geometries can use
 * the empty default implementation.
 */
inline void RigidBody::handleRotation()
{
   WALBERLA_ABORT( "Invalid call of default handle function" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Handling a fixation change of a contained subordinate body.
 *
 * \return void
 *
 * This function handles a fixation change of one of the contained subordinate bodies. Derived
 * compound geometries that contain other primitive bodies are required to override this function
 * in order to react to the fixation change. All primitive geometries can use the empty default
 * implementation.
 */
inline void RigidBody::handleFixation()
{
   WALBERLA_ABORT( "Invalid call of default handle function" );
}
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
