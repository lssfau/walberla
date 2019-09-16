#//======================================================================================================================
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
//! \file HCSITS.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the hard contact solver
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/cr/ICR.h"
#include "pe/communication/PackNotification.h"
#include "pe/communication/RigidBodyVelocityCorrectionNotification.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/Types.h"

#include <core/config/Config.h>
#include <core/DataTypes.h>
#include <core/debug/Debug.h>
#include <core/logging/Logging.h>
#include <core/math/Matrix2.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector2.h>
#include <core/math/Vector3.h>
#include <core/mpi/BufferSystem.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>
#include <core/NonCopyable.h>

namespace walberla {
namespace pe {
namespace cr {


//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Particular implementation of the collision resoution for the hard contacts.
 *
 * The following code example illustrates the setup of the solver:
 * \snippet PeDocumentationSnippets.cpp Setup HCSITS
 */
class HardContactSemiImplicitTimesteppingSolvers
   : public ICR
   , private NonCopyable
{
private:
   // bodies
   struct BodyCache
   {
      void resize(const size_t n);

      std::vector<Vec3> v_, w_, dv_, dw_;
   };
   std::map<IBlockID::IDType, BodyCache> blockToBodyCache_;

   // contacts
   struct ContactCache
   {
      void resize(const size_t n);

      std::vector<Vec3>   r1_, r2_; // vector pointing from body1/body2 to the contact position
      std::vector<BodyID> body1_, body2_;
      std::vector<Vec3>   n_, t_, o_; // contact normal and the two other directions perpendicular to the normal
      std::vector<real_t> dist_; // overlap length, a contact is present if dist_ < 0
      std::vector<real_t> mu_; // contact friction
      std::vector<Mat3>   diag_nto_; 
      std::vector<Mat3>   diag_nto_inv_;
      std::vector<Mat2>   diag_to_inv_;
      std::vector<real_t> diag_n_inv_;
      std::vector<Vec3>   p_;
   };
   std::map<IBlockID::IDType, ContactCache> blockToContactCache_;

public:
   //**Definition of relaxation models ************************************************************
   enum RelaxationModel {
      InelasticFrictionlessContact,
      ApproximateInelasticCoulombContactByDecoupling,
      ApproximateInelasticCoulombContactByOrthogonalProjections,
      InelasticCoulombContactByDecoupling,
      InelasticCoulombContactByOrthogonalProjections,
      InelasticGeneralizedMaximumDissipationContact
   };
   //**********************************************************************************************
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit HardContactSemiImplicitTimesteppingSolvers( const shared_ptr<BodyStorage>&      globalBodyStorage,
                                                        const shared_ptr<BlockStorage>&     blockStorage,
                                                        domain_decomposition::BlockDataID   storageID,
                                                        domain_decomposition::BlockDataID   ccdID,
                                                        domain_decomposition::BlockDataID   fcdID,
                                                        WcTimingTree* tt = NULL );
   //@}
   //**********************************************************************************************
   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~HardContactSemiImplicitTimesteppingSolvers();
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   virtual inline real_t            getMaximumPenetration()        const WALBERLA_OVERRIDE;
   virtual inline size_t            getNumberOfContacts()          const WALBERLA_OVERRIDE;
   virtual inline size_t            getNumberOfContactsTreated()   const WALBERLA_OVERRIDE;
   inline const std::map<IBlockID::IDType, ContactCache> getContactCache() const { return blockToContactCache_; }
   inline real_t                    getSpeedLimitFactor() const;
   inline size_t                    getMaxIterations() const { return maxIterations_; }
   inline real_t                    getOverRelaxationParameter() const { return overRelaxationParam_; }
   inline real_t                    getRelaxationParameter() const { return relaxationParam_; }
   inline real_t                    getErrorReductionParameter() const { return erp_; }
   inline RelaxationModel           getRelaxationModel() const { return relaxationModel_; }
   //@}
   //**********************************************************************************************

   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   inline void            setOverRelaxationParameter( real_t omega );
   inline void            setRelaxationParameter( real_t f );
   inline void            setMaxIterations( size_t n );
   inline void            setRelaxationModel( RelaxationModel relaxationModel );
   inline void            setErrorReductionParameter( real_t erp );
   inline void            setAbortThreshold( real_t threshold );
   inline void            setSpeedLimiter( bool active, const real_t speedLimitFactor = real_t(0.0) );
   //@}
   //**********************************************************************************************

   //**Query functions*****************************************************************************
   /*!\name Query functions */
   //@{
   inline bool            isSyncRequired()        const;
   inline bool            isSyncRequiredLocally() const;
   inline bool            isSpeedLimiterActive() const;
   //@}
   //**********************************************************************************************

   /// forwards to timestep
   /// Convenience operator to make class a functor.
   void operator()(const real_t dt) { timestep(dt); }
   /// Advances the simulation dt seconds.
   void timestep( const real_t dt );

private:

   //**Simulation functions************************************************************************
   /*!\name Simulation functions */
   //@{
   void resolveContacts( const Contacts& contacts, real_t dt );
   real_t relaxInelasticFrictionlessContacts( real_t dtinv,
                                              HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                              HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache );
   real_t relaxApproximateInelasticCoulombContactsByDecoupling( real_t dtinv,
                                                                HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache );
   real_t relaxInelasticCoulombContactsByDecoupling( real_t dtinv,
                                                     HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                     HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache );
   real_t relaxInelasticCoulombContactsByOrthogonalProjections( real_t dtinv,
                                                                bool approximate,
                                                                HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache );
   real_t relaxInelasticGeneralizedMaximumDissipationContacts( real_t dtinv,
                                                               HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                               HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache );
   //@}
   //**********************************************************************************************

   //**Communication functions*********************************************************************
   /*!\name Communication functions */
   //@{
   void parseVelocityCorrection(mpi::RecvBuffer& rb, BodyStorage& bodyStorage, BodyCache& bodyCache);
   void parseVelocityCorrectionShadow(mpi::RecvBuffer& rb, BodyStorage& bodyStorage, BodyCache& bodyCache);
   void synchronizeVelocities();
   //@}
   //**********************************************************************************************

   //**Time-integration functions******************************************************************
   /*!\name Time-integration functions */
   //@{
   void initializeVelocityCorrections( BodyID body, Vec3& dv, Vec3& dw, real_t dt ) const;
   void integratePositions( BodyID body, Vec3 v, Vec3 w, real_t dt ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
          bool checkUpdateFlags();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   shared_ptr<BodyStorage>           globalBodyStorage_;
   shared_ptr<BlockStorage>          blockStorage_;
   domain_decomposition::BlockDataID storageID_;
   domain_decomposition::BlockDataID ccdID_;
   domain_decomposition::BlockDataID fcdID_;
   WcTimingTree*                     tt_;


   real_t erp_;                       //!< The error reduction parameter (0 <= erp_ <= 1).
   size_t maxIterations_;             //!< Maximum number of iterations.
   size_t iteration_;
   size_t maxSubIterations_;          //!< Maximum number of iterations of iterative solvers in the one-contact problem.
   real_t abortThreshold_;            //!< If L-infinity iterate difference drops below this threshold the iteration is aborted.
   RelaxationModel relaxationModel_;  //!< The method used to relax unilateral contacts
   real_t overRelaxationParam_;       //!< Parameter specifying the convergence speed for othogonal projection models.
   real_t relaxationParam_;           //!< Parameter specifying underrelaxation of velocity corrections for boundary bodies.
   real_t maximumPenetration_;
   size_t numContacts_;
   size_t numContactsTreated_;

   bool   speedLimiterActive_;        //!< is the speed limiter active?
   real_t speedLimitFactor_;          //!< what multiple of boundingbox edge length is the body allowed to travel in one timestep

   //**********************************************************************************************
   /*! \cond WALBERLA_INTERNAL */
   /*!\brief Functor for comparing the system ID of two bodies.
    *
    * Returns true if the system ID of the first body is less than the system ID of the second body.
    */
   struct LessSystemID {
      //**Binary function call operator************************************************************
      /*!\name Binary function call operator */
      //@{
      bool operator()( BodyID b1, BodyID b2 ) const {
         return b1->getSystemID() < b2->getSystemID();
      }
      //@}
      //*******************************************************************************************
   };
   /*! \endcond */
   //**********************************************************************************************
   bool requireSync_;         //!< Flag indicating whether this process requires a synchronization prior to the next time step.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************

typedef HardContactSemiImplicitTimesteppingSolvers HCSITS;

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Returns the maximum penetration depth found in the last collision detection.
 *
 * \return The maximum penetration depth found in the last collision detection.
 *
 * Only contacts treated on the local process are considered.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::getMaximumPenetration() const
{
   return maximumPenetration_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of contacts found by the last collision detection.
 *
 * \return The number of contacts found by the last collision detection.
 *
 * Only contacts treated on the local process are counted.
 */
inline size_t HardContactSemiImplicitTimesteppingSolvers::getNumberOfContacts() const
{
   return numContacts_;
}

inline size_t HardContactSemiImplicitTimesteppingSolvers::getNumberOfContactsTreated() const
{
   return numContactsTreated_;
}

inline real_t HardContactSemiImplicitTimesteppingSolvers::getSpeedLimitFactor() const
{
   return speedLimitFactor_;
}
//*************************************************************************************************




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Sets the relaxation parameter for boundary bodies.
 *
 * \param f The overrelaxation parameter.
 * \return void
 *
 * The overrelaxation parameter \omega is only used when the relaxation model is one of
 * - ApproximateInelasticCoulombContactByOrthogonalProjections
 * - InelasticCoulombContactByOrthogonalProjections
 *
 * It is used to control the convergence of the model. Large values show faster convergence,
 * but they can also lead to divergence ("exploding" particles). The default values is 1.0.
 */
inline void HardContactSemiImplicitTimesteppingSolvers::setOverRelaxationParameter( real_t omega )
{
   WALBERLA_ASSERT_GREATER( omega, 0, "Overrelaxation parameter must be positive." );

   overRelaxationParam_ = omega;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets the relaxation parameter for boundary bodies.
 *
 * \param f The relaxation parameter.
 * \return void
 *
 * The iterative solvers are a mixture of non-linear Gauss-Seidel and Jacobi solvers. This might
 * require underrelaxation. The parameter must be positive. Note that for dilute systems the
 * solver might need stronger underrelaxation (smaller \a f) than for dense systems.
 */
inline void HardContactSemiImplicitTimesteppingSolvers::setRelaxationParameter( real_t f )
{
   WALBERLA_ASSERT_GREATER( f, 0, "Relaxation parameter must be positive." );

   relaxationParam_ = f;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets the maximum number of iterations performed by the iterative solver.
 *
 * \param n The maximum  number of iterations.
 * \return void
 */
inline void HardContactSemiImplicitTimesteppingSolvers::setMaxIterations( size_t n )
{
   maxIterations_ = n;
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Sets the relaxation model used by the iterative solver.
 *
 * \param relaxationModel The relaxation model to be used by the iterative solver.
 * \return void
 */
inline void HardContactSemiImplicitTimesteppingSolvers::setRelaxationModel( RelaxationModel relaxationModel )
{
   relaxationModel_ = relaxationModel;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets the error reduction parameter.
 *
 * \param erp The error reduction parameter.
 * \return void
 *
 * If body shapes overlap by x at a contact then the contact resolution aims to remove erp*x of the
 * overlap. Thus the error reduction parameter must be between 0 and 1. 0 corresponds to no error
 * reduction and is the default. 1 corresponds to full error reduction. Note that error reduction
 * (constraint stabilization) introduces additional energy to the system.
 */
inline void HardContactSemiImplicitTimesteppingSolvers::setErrorReductionParameter( real_t erp )
{
   WALBERLA_ASSERT_GREATER_EQUAL( erp, 0, "Error reduction parameter out of range." );
   WALBERLA_ASSERT_LESS_EQUAL( erp, 1, "Error reduction parameter out of range." );

   erp_ = erp;
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Sets the threshold of movement of a particle during collsion resolulution
*
* \param threshold If movement is smaller than threshold, col. resolution is stopped.
* \return void
*/
inline void HardContactSemiImplicitTimesteppingSolvers::setAbortThreshold( real_t threshold )
{
   abortThreshold_ = threshold;
}
//*************************************************************************************************




//*************************************************************************************************
/*!\brief Activates/Deactivates the speed limiter and sets the limit
*
* \param active activate/deactivate speed limtier
* \param speedLimitFactor size of bounding box will be multiplied by this factor to get the maximal distance a body is allowed to travel within one timestep
* \return void
*/
inline void HardContactSemiImplicitTimesteppingSolvers::setSpeedLimiter( bool active, const real_t speedLimitFactor )
{
   speedLimiterActive_ = active;
   speedLimitFactor_   = speedLimitFactor;
}
//*************************************************************************************************




//=================================================================================================
//
//  QUERY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns if a synchronization is required before the next time step.
 *
 * \return True if a synchronization is strictly required, false if a synchronization is possibly
 *         required for correct dynamics but not enforced.
 *
 * Insertion or removal of bodies from the simulation \em can require a subsequent synchronization
 * before performing the next time step. If e.g. no other process obtained a shadow copy of a
 * body to be removed then a synchronization is \em not enforced. However, if a neighbor has a
 * shadow copy a synchronization is required. Changing e.g. velocities or positions can lead to
 * inconsistent descriptions of bodies across process boundaries but synchronization is not
 * enforced. In this case it is the users obligation to synchronize whenever necessary.
 *
 * WARNING: This query function uses an expensive allreduce MPI operation to determine the
 * result!
 */
inline bool HardContactSemiImplicitTimesteppingSolvers::isSyncRequired() const
{
   WALBERLA_MPI_SECTION()
   {
      // No synchronization ever necessary if we compute on a single process
      if( mpi::MPIManager::instance()->numProcesses() <= 1 )
         return false;

      char requireSync( requireSync_ );

      mpi::allReduceInplace(requireSync, mpi::LOGICAL_OR);
      return requireSync != 0;
   } else
   {
      return false;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns if a synchronization is required by the local process before the next time step.
 *
 * \return True if a synchronization is strictly required by the local process, false otherwise.
 *
 * Insertion or removal of bodies from the simulation \em can require a subsequent synchronization
 * before performing the next time step. If e.g. no other process obtained a shadow copy of a
 * body to be removed then a synchronization is \em not enforced. However, if a neighbor has a
 * shadow copy a synchronization is required. Changing e.g. velocities or positions can lead to
 * inconsistent descriptions of bodies across process boundaries but synchronization is not
 * enforced. In this case it is the users obligation to synchronize whenever necessary.
 */
inline bool HardContactSemiImplicitTimesteppingSolvers::isSyncRequiredLocally() const
{
   return requireSync_;
}
//*************************************************************************************************




//*************************************************************************************************
/*!\brief Returns if speed limiter is currently active and working.
 *
 * \return status of the speed limiter
 */
inline bool HardContactSemiImplicitTimesteppingSolvers::isSpeedLimiterActive() const
{
   return speedLimiterActive_;
}
//*************************************************************************************************


} // namespace cr
} // namespace pe

/**
 * \brief configures HardContactSemiImplicitTimesteppingSolvers with parameters from config file
 * \param config handle to config block
 * \param cr collision resolution object to configure
 */
void configure( const Config::BlockHandle& config, pe::cr::HCSITS& cr);

} // namespace walberla

#include "HCSITS.impl.h"
