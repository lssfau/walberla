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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/IContactAccessor.h>
#include <mesa_pd/data/Flags.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/math/Limits.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
/**
 * Apply a HCSITS (Hard-contact semi-implicit timestepping solver) relaxation to all contacts.
 * Call this kernel on all bodies to perform on relaxation step of the solver.
 * There exist different friction and relaxation models. See comments on their respective methods
 * for an extended documentation. Use setRelaxationModel() to set it.
 * getDeltaMax() can be used to query the maximum change in the impulses applied wrt to the last iteration.
 *
 * Use setMaxSubIterations() to set the maximum number of iterations performed to solve optimization problems arising
 * during the relaxation (in InelasticGeneralizedMaximumDissipationContact and InelasticCoulombContactByDecoupling).
 *
 * Use setCor() to set the coefficient of restitution between 0 and 1 for the InelasticProjectedGaussSeidel model.
 * All other models describe inelastic contacts.
 *
 * Example of Using the HCSITS kernels in a simulation:
 * \code
 *     // ps: Particle storage
 *     // cs: Contact storage
 *     // pa: Particle accessor
 *     // ca: Contact accessor
 *
 *     // Detect contacts first and fill contact storage
 *
 *     mesa_pd::mpi::ReduceProperty reductionKernel;
 *     mesa_pd::mpi::BroadcastProperty broadcastKernel;
 *
 *     kernel::IntegrateBodiesHCSITS integration;
 *
 *     kernel::InitContactsForHCSITS initContacts(1);
 *     initContacts.setFriction(0,0,real_t(0.2));
 *     kernel::InitBodiesForHCSITS initBodies;
 *
 *     kernel::HCSITSRelaxationStep relaxationStep;
 *     relaxationStep.setCor(real_t(0.1)); // Only effective for PGSM
 *
 *     // Init Contacts and Bodies (order is arbitrary)
 *     cs.forEachContact(false, kernel::SelectAll(), ca, initContacts, ca, pa);
 *     ps.forEachParticle(false, kernel::SelectAll(), pa, initBodies, pa, dt);
 *
 *     // Reduce and Broadcast velocities with relaxation parameter 1 before the iteration
 *     VelocityUpdateNotification::Parameters::relaxationParam = real_t(1.0);
 *     reductionKernel.operator()<VelocityCorrectionNotification>(*ps);
 *     broadcastKernel.operator()<VelocityUpdateNotification>(*ps);
 *
 *     VelocityUpdateNotification::Parameters::relaxationParam = real_t(0.8);
 *     for(int j = 0; j < 10; j++){
 *        cs.forEachContact(false, kernel::SelectAll(), ca, relaxationStep, ca, pa, dt);
 *        reductionKernel.operator()<VelocityCorrectionNotification>(*ps);
 *        broadcastKernel.operator()<VelocityUpdateNotification>(*ps);
 *     }
 *
 *     ps.forEachParticle(false, kernel::SelectAll(), pa, integration, pa, dt);
 * \endcode
 * \ingroup mesa_pd_kernel
 *
 *
 */
class HCSITSRelaxationStep
{
public:
   enum RelaxationModel {
      InelasticFrictionlessContact,
      ApproximateInelasticCoulombContactByDecoupling,
      ApproximateInelasticCoulombContactByOrthogonalProjections,
      InelasticCoulombContactByDecoupling,
      InelasticCoulombContactByOrthogonalProjections,
      InelasticGeneralizedMaximumDissipationContact,
      InelasticProjectedGaussSeidel
   };
   
   // Default constructor sets the default values
   HCSITSRelaxationStep() :
   {%- for prop in properties %}
   {{prop.name}}_( {{prop.defValue}} ){{ "," if not loop.last}}
   {%- endfor %}
   {}

   // Getter and Setter Functions
   {%- for prop in properties %}
   inline const {{prop.type}}& get{{prop.name | capFirst}}() const {return {{prop.name}}_;}
   inline void set{{prop.name | capFirst}}({{prop.type}} v) { {{prop.name}}_ = v;}
   {% endfor %}

   /**
    * Perform relaxation for a single contact.
    * \param cid The index of the contact with the accessor ca.
    * \param ca The contact accessor.
    * \param pa The particle accessor.
    * \param dt The timestep size used.
    */
   template <typename CAccessor, typename PAccessor>
   void operator()(size_t cid, CAccessor &ca, PAccessor& pa,  real_t dt);

   template <typename CAccessor, typename PAccessor>
   real_t relaxInelasticFrictionlessContacts(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa);

   template <typename CAccessor, typename PAccessor>
   real_t relaxApproximateInelasticCoulombContactsByDecoupling(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa);

   template <typename CAccessor, typename PAccessor>
   real_t relaxInelasticCoulombContactsByDecoupling(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa );

   template <typename CAccessor, typename PAccessor>
   real_t relaxInelasticCoulombContactsByOrthogonalProjections(size_t cid, real_t dtinv, bool approximate, CAccessor &ca, PAccessor& pa );

   template <typename CAccessor, typename PAccessor>
   real_t relaxInelasticGeneralizedMaximumDissipationContacts(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa );

   template <typename CAccessor, typename PAccessor>
   real_t relaxInelasticContactsByProjectedGaussSeidel(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa );

private:
   // List properties
   {% for prop in properties %}
   {{prop.type}} {{prop.name }}_;
   {%- endfor %}
};


template <typename CAccessor, typename PAccessor>
inline void HCSITSRelaxationStep::operator()(size_t cid, CAccessor &ca, PAccessor& pa,  real_t dt)
{
   static_assert(std::is_base_of<data::IContactAccessor, CAccessor>::value, "please provide a valid contact accessor");
   static_assert(std::is_base_of<data::IAccessor, PAccessor>::value, "please provide a valid accessor");

   real_t dtinv(real_t(1.0)/dt);
   switch( getRelaxationModel() ) {
      case InelasticFrictionlessContact:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticFrictionlessContacts(cid, dtinv, ca, pa )));
         break;

      case ApproximateInelasticCoulombContactByDecoupling:
         setDeltaMax(std::max( getDeltaMax(), relaxApproximateInelasticCoulombContactsByDecoupling( cid, dtinv, ca, pa )));
         break;

      case ApproximateInelasticCoulombContactByOrthogonalProjections:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticCoulombContactsByOrthogonalProjections( cid, dtinv, true, ca, pa )));
         break;

      case InelasticCoulombContactByDecoupling:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticCoulombContactsByDecoupling( cid, dtinv, ca, pa )));
         break;

      case InelasticCoulombContactByOrthogonalProjections:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticCoulombContactsByOrthogonalProjections( cid, dtinv, false, ca, pa )));
         break;

      case InelasticGeneralizedMaximumDissipationContact:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticGeneralizedMaximumDissipationContacts( cid, dtinv, ca, pa )));
         break;

      case InelasticProjectedGaussSeidel:
         setDeltaMax(std::max( getDeltaMax(), relaxInelasticContactsByProjectedGaussSeidel( cid, dtinv, ca, pa )));
         break;

      default:
         throw std::runtime_error( "Unsupported relaxation model." );
}

WALBERLA_LOG_DETAIL("Delta Max: " << getDeltaMax());

}

//*************************************************************************************************
/*!\brief Relaxes The contact with ID cid once. The contact model is for inelastic unilateral contacts without friction.
 *
 * \param cid The index of the contact
 * \param ca The contact accessor
 * \param pa The particle accessor
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * persisting solutions if valid.
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxInelasticFrictionlessContacts(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa)
{
   real_t delta_max( 0 );

   // Remove velocity corrections of this contact's reaction.
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));

   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   Vec3 gdot    ( ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) -
                  ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) +
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) -
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) /* + diag_[i] * p */ );


   // Change from the global world frame to the contact frame
   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
   Vec3 gdot_nto( contactframe.getTranspose() * gdot );

   // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
   gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + ca.getDistance(cid) ) * dtinv;

   if( gdot_nto[0] >= 0 ) {
      // Contact is separating if no contact reaction is present at contact i.

      delta_max = std::max( delta_max, std::max( std::abs( ca.getP(cid)[0] ), std::max( std::abs( ca.getP(cid)[1] ), std::abs( ca.getP(cid)[2] ) ) ) );
      ca.getPRef(cid) = Vec3();

      // No need to apply zero impulse.
   }
   else {
      // Contact is persisting.

      // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
      Vec3 p_wf( ca.getNormal(cid) * ( -ca.getDiag_n_inv(cid) * gdot_nto[0] ) );
      Vec3 dp( ca.getP(cid) - p_wf );
      delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

      ca.getPRef(cid) = p_wf;

      // Apply impulse right away.
      pa.getDvRef(bId1) += pa.getInvMass(bId1) * ca.getP(cid);
      pa.getDwRef(bId1) += pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
      pa.getDvRef(bId2) -= pa.getInvMass(bId2) * ca.getP(cid);
      pa.getDwRef(bId2) -= pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));
   }

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with approximate Coulomb friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * other solutions if valid. Static solutions are preferred over dynamic solutions. Dynamic
 * solutions are computed by decoupling the normal from the frictional components. That is
 * for a dynamic contact the normal component is relaxed first followed by the frictional
 * components. The determination of the frictional components does not perform any subiterations
 * and guarantees that the friction partially opposes slip.
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxApproximateInelasticCoulombContactsByDecoupling(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa )
{
   real_t delta_max( 0 );

   // Remove velocity corrections of this contact's reaction.
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));


   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   Vec3 gdot    ( ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) -
                  ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) +
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) -
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) /* + diag_[i] * p */ );

   // Change from the global world frame to the contact frame
   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
   Vec3 gdot_nto( contactframe.getTranspose() * gdot );

   // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
   gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + ca.getDistance(cid) ) * dtinv;


   if( gdot_nto[0] >= 0 ) {
      // Contact is separating if no contact reaction is present at contact i.
      delta_max = std::max( delta_max, std::max( std::abs( ca.getP(cid)[0] ), std::max( std::abs( ca.getP(cid)[1] ), std::abs( ca.getP(cid)[2] ) ) ) );
      ca.getPRef(cid) = Vec3();

      // No need to apply zero impulse.
   }
   else {
      // Contact is persisting (either static or dynamic).

      // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
      //Vec3 p_cf( -( contactCache.diag_nto_inv_[i] * gdot_nto ) );
      Vec3 p_cf( -( ca.getDiag_nto_inv(cid) * gdot_nto ) );

      // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
      // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
      // b = [0.01 -1 -1]';
      // A\b    \approx [-0.19 -2.28 -1.21]'
      // eig(A) \approx [ 0.40  0.56  1.04]'

      //real_t flimit( contactCache.mu_[i] * p_cf[0] );
      real_t flimit( ca.getMu(cid) * p_cf[0] );
      real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
      if( fsq > flimit * flimit || p_cf[0] < 0 ) {
         // Contact cannot be static so it must be dynamic.
         // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

         // For simplicity we change to a simpler relaxation scheme here:
         // 1. Relax normal reaction with the tangential components equal to the previous values
         // 2. Relax tangential components with the newly relaxed normal reaction
         // Note: The better approach would be to solve the true 3x3 block problem!
         // Warning: Simply projecting the frictional components is wrong since then the normal action is no longer 0 and simulations break.

         // Add the action of the frictional reactions from the last iteration to the relative contact velocity in normal direction so we can relax it separately.
         // TODO This can be simplified:
         //p_cf = trans( contactframe ) * contactCache.p_[i];
         //p_cf[0] = 0;
         //p_[i] = contactframe * p_cf;


         //Vec3 p_tmp = ( contactCache.t_[i] * contactCache.p_[i] ) * contactCache.t_[i] + ( contactCache.o_[i] * contactCache.p_[i] ) * contactCache.o_[i];
         Vec3 p_tmp = ( ca.getT(cid) * ca.getP(cid)  ) * ca.getT(cid)  + ( ca.getO(cid) * ca.getP(cid) ) * ca.getO(cid);


         //real_t gdot_n = gdot_nto[0] + contactCache.n_[i] * ( ( contactCache.body1_[i]->getInvInertiaBF() * ( contactCache.r1_[i] % p_tmp ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertiaBF() * ( contactCache.r2_[i] % p_tmp ) ) % contactCache.r2_[i] /* + diag_[i] * p */ );
         real_t gdot_n = gdot_nto[0] + ca.getNormal(cid) * ( ( pa.getInvInertiaBF(bId1) * ( ca.getR1(cid) % p_tmp ) ) % ca.getR1(cid) + ( pa.getInvInertiaBF(bId2) * ( ca.getR2(cid) % p_tmp ) ) % ca.getR2(cid) /* + diag_[i] * p */ );

         //p_cf[0] = -( contactCache.diag_n_inv_[i] * gdot_n );
         p_cf[0] = -( ca.getDiag_n_inv(cid) * gdot_n );

         // We cannot be sure that gdot_n <= 0 here and thus p_cf[0] >= 0 since we just modified it with the old values of the tangential reactions! => Project
         p_cf[0] = std::max( real_c( 0 ), p_cf[0] );

         // Now add the action of the normal reaction to the relative contact velocity in the tangential directions so we can relax the frictional components separately.
         p_tmp = ca.getNormal(cid)  * p_cf[0];
         Vec3 gdot2 = gdot + ( pa.getInvInertiaBF(bId1) * ( ca.getR1(cid) % p_tmp ) ) % ca.getR1(cid) + ( pa.getInvInertiaBF(bId2)* ( ca.getR2(cid) % p_tmp ) ) % ca.getR2(cid);
         Vec2 gdot_to;
         gdot_to[0] = ca.getT(cid) * gdot2;
         gdot_to[1] = ca.getO(cid) * gdot2;

         //Vec2 ret = -( contactCache.diag_to_inv_[i] * gdot_to );
         Vec2 ret = -( ca.getDiag_to_inv(cid) * gdot_to );
         p_cf[1] = ret[0];
         p_cf[2] = ret[1];

         flimit = ca.getMu(cid) * p_cf[0];
         fsq = p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2];
         if( fsq > flimit * flimit ) {
            const real_t f( flimit / std::sqrt( fsq ) );
            p_cf[1] *= f;
            p_cf[2] *= f;
         }
      }
      else {
         // Contact is static.
      }
      Vec3 p_wf( contactframe * p_cf );

      Vec3 dp( ca.getP(cid) - p_wf );
      delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

      ca.getPRef(cid) = p_wf;

      // Apply impulse right away.
      pa.getDvRef(bId1) += pa.getInvMass(bId1) * ca.getP(cid);
      pa.getDwRef(bId1) += pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
      pa.getDvRef(bId2) -= pa.getInvMass(bId2) * ca.getP(cid);
      pa.getDwRef(bId2) -= pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));
   }
   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************

/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with Coulomb friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * other solutions if valid. Static solutions are preferred over dynamic solutions. Dynamic
 * solutions are computed by decoupling the normal from the frictional components. That is
 * for a dynamic contact the normal component is relaxed first followed by the frictional
 * components. How much the frictional components directly oppose slip as required by the Coulomb
 * friction model depends on the number of subiterations performed. If no subiterations are
 * performed the friction is guaranteed to be at least partially dissipative.
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxInelasticCoulombContactsByDecoupling(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa )
{
   real_t delta_max( 0 );

   // Remove velocity corrections of this contact's reaction.
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));


   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   Vec3 gdot    ( ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) -
                  ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) +
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) -
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) /* + diag_[i] * p */ );

   // Change from the global world frame to the contact frame
   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
   Vec3 gdot_nto( contactframe.getTranspose() * gdot );

   // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
   gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + ca.getDistance(cid) ) * dtinv;

   //WALBERLA_LOG_WARNING( "Contact #" << i << " is\nA = \n" << contactCache.diag_nto_[i] << "\nb = \n" << gdot_nto << "\nmu = " << contactCache.mu_[i] );

   if( gdot_nto[0] >= 0 ) {
      // Contact is separating if no contact reaction is present at contact i.

      delta_max = std::max( delta_max, std::max( std::abs( ca.getP(cid)[0] ), std::max( std::abs( ca.getP(cid)[1] ), std::abs( ca.getP(cid)[2] ) ) ) );
      ca.getPRef(cid) = Vec3();
      //WALBERLA_LOG_WARNING( "Contact #" << i << " is separating." );

      // No need to apply zero impulse.
   }
   else {
      // Contact is persisting (either static or dynamic).

      // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
      Vec3 p_cf( -( ca.getDiag_nto_inv(cid) * gdot_nto ) );

      // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
      // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
      // b = [0.01 -1 -1]';
      // A\b    \approx [-0.19 -2.28 -1.21]'
      // eig(A) \approx [ 0.40  0.56  1.04]'

      real_t flimit( ca.getMu(cid) * p_cf[0] );
      real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
      if( fsq > flimit * flimit || p_cf[0] < 0 ) {
         // Contact cannot be static so it must be dynamic.
         // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

         for (int j = 0; j < 20; ++j) {
            // For simplicity we change to a simpler relaxation scheme here:
            // 1. Relax normal reaction with the tangential components equal to the previous values
            // 2. Relax tangential components with the newly relaxed normal reaction
            // Note: The better approach would be to solve the true 3x3 block problem!
            // Warning: Simply projecting the frictional components is wrong since then the normal action is no longer 0 and simulations break.

            Vec3 gdotCorrected;
            real_t gdotCorrected_n;
            Vec2 gdotCorrected_to;

            // Calculate the relative contact velocity in the global world frame (if no normal contact reaction is present at contact i)
            p_cf[0] = 0;
            //                       |<- p_cf is orthogonal to the normal and drops out in next line ->|
            gdotCorrected   = /* ( contactCache.body1_[i]->getInvMass() + contactCache.body2_[i]->getInvMass() ) * p_cf  */ gdot + ( pa.getInvInertiaBF(bId1) * ( ca.getR1(cid) % ( ca.getT(cid) * p_cf[1] + ca.getO(cid) * p_cf[2] ) ) ) % ca.getR1(cid) + ( pa.getInvInertiaBF(bId2) * ( ca.getR2(cid) % ( ca.getT(cid) * p_cf[1] + ca.getO(cid) * p_cf[2] ) ) ) % ca.getR2(cid);
            gdotCorrected_n = ca.getNormal(cid) * gdotCorrected + ca.getDistance(cid) * dtinv;

            // Relax normal component.
            p_cf[0] = std::max( real_c( 0 ), -( ca.getDiag_n_inv(cid) * gdotCorrected_n ) );

            // Calculate the relative contact velocity in the global world frame (if no frictional contact reaction is present at contact i)
            p_cf[1] = p_cf[2] = real_c( 0 );
            //                       |<- p_cf is orthogonal to the tangential plane and drops out   ->|
            gdotCorrected   = /* ( contactCache.body1_[i]->getInvMass() + contactCache.body2_[i]->getInvMass() ) * p_cf */ gdot + ( pa.getInvInertiaBF(bId1) * ( ca.getR1(cid) % ( ca.getNormal(cid) * p_cf[0] ) ) ) % ca.getR1(cid) + ( pa.getInvInertiaBF(bId2) * ( ca.getR2(cid) % ( ca.getNormal(cid) * p_cf[0] ) ) ) % ca.getR2(cid);
            gdotCorrected_to[0] = ca.getT(cid) * gdotCorrected;
            gdotCorrected_to[1] = ca.getO(cid) * gdotCorrected;

            // Relax frictional components.
            Vec2 ret = -( ca.getDiag_to_inv(cid) * gdotCorrected_to );
            p_cf[1] = ret[0];
            p_cf[2] = ret[1];

            flimit = ca.getMu(cid) * p_cf[0];
            fsq = p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2];
            if( fsq > flimit * flimit ) {
               // 3.2.1 Decoupling
               // \tilde{x}^0 = p_cf[1..2]

               // Determine \tilde{A}
               Mat2 diag_to( ca.getDiag_nto(cid)(1, 1), ca.getDiag_nto(cid)(1, 2), ca.getDiag_nto(cid)(2, 1), ca.getDiag_nto(cid)(2, 2) );

               const real_t f( flimit / std::sqrt( fsq ) );
               //p_cf[1] *= f;
               //p_cf[2] *= f;

               // Determine search interval for Golden Section Search
               const real_t phi( real_c(0.5) * ( real_c(1) + std::sqrt( real_c( 5 ) ) ) );
               real_t shift( std::atan2( -p_cf[2], p_cf[1] ) );
               real_t acos_f( std::acos( f ) );

               //WALBERLA_LOG_WARNING( acos_f << " " << shift );

               real_t alpha_left( -acos_f - shift );
               //Vec2 x_left( flimit * std::cos( alpha_left ), flimit * std::sin( alpha_left ) );
               //real_t f_left( 0.5 * trans( x_left ) * ( diag_to * x_left ) - trans( x_left ) * ( -gdot_to ) );

               real_t alpha_right( acos_f - shift );
               //Vec2 x_right( flimit * std::cos( alpha_right ), flimit * std::sin( alpha_right ) );
               //real_t f_right( 0.5 * trans( x_right ) * ( diag_to * x_right ) - trans( x_right ) * ( -gdot_to ) );

               real_t alpha_mid( ( alpha_right + alpha_left * phi ) / ( 1 + phi ) );
               Vec2 x_mid( flimit * std::cos( alpha_mid ), flimit * std::sin( alpha_mid ) );
               real_t f_mid( real_c(0.5) * x_mid * ( diag_to * x_mid ) - x_mid * ( -gdotCorrected_to ) );

               bool leftlarger = false;
               for( size_t k = 0; k < getMaxSubIterations(); ++k ) {
                  real_t alpha_next( alpha_left + ( alpha_right - alpha_mid ) );
                  Vec2 x_next( flimit * std::cos( alpha_next ), flimit * std::sin( alpha_next ) );
                  real_t f_next( real_c(0.5) * x_next * ( diag_to * x_next ) - x_next * ( -gdotCorrected_to ) );
                  //WALBERLA_LOG_WARNING( "[(" << alpha_left << ", ?); (" << alpha_mid << ", " << f_mid << "); (" << alpha_right << ", ?)] <- (" << alpha_next << ", " << f_next << ")" );
                  //WALBERLA_LOG_WARNING( "left: " << alpha_mid - alpha_left << "  right: " << alpha_right - alpha_mid << "  ll: " << leftlarger );
                  //WALBERLA_ASSERT(leftlarger ? (alpha_mid - alpha_left > alpha_right - alpha_mid) : (alpha_mid - alpha_left < alpha_right - alpha_mid), "ll inconsistent!" );

                  if (leftlarger) {
                     // left interval larger
                     if( f_next < f_mid ) {
                        alpha_right = alpha_mid;
                        alpha_mid   = alpha_next;
                        x_mid       = x_next;
                        f_mid       = f_next;
                        leftlarger = true;
                     }
                     else {
                        alpha_left  = alpha_next;
                        leftlarger = false;
                     }
                  }
                  else {
                     // right interval larger
                     if( f_next < f_mid ) {
                        alpha_left = alpha_mid;
                        alpha_mid  = alpha_next;
                        x_mid      = x_next;
                        f_mid      = f_next;
                        leftlarger = false;
                     }
                     else {
                        alpha_right = alpha_next;
                        leftlarger = true;
                     }
                  }
               }
               //WALBERLA_LOG_WARNING( "dalpha = " << alpha_right - alpha_left );

               p_cf[1] = x_mid[0];
               p_cf[2] = x_mid[1];
            }
         }
         //WALBERLA_LOG_WARNING( "Contact #" << i << " is dynamic." );
      }
      else {
         // Contact is static.
         //WALBERLA_LOG_WARNING( "Contact #" << i << " is static." );
      }

      //WALBERLA_LOG_WARNING( "Contact reaction in contact frame: " << p_cf << "\n" << ca.getDiag_nto(cid)*p_cf + gdot_nto );
      Vec3 p_wf( contactframe * p_cf );
      Vec3 dp( ca.getP(cid) - p_wf );
      delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

      ca.getPRef(cid) = p_wf;

      // Apply impulse right away.
      pa.getDvRef(bId1) += pa.getInvMass(bId1) * ca.getP(cid);
      pa.getDwRef(bId1) += pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
      pa.getDvRef(bId2) -= pa.getInvMass(bId2) * ca.getP(cid);
      pa.getDwRef(bId2) -= pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));
   }

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with Coulomb friction.
 *
 * \param dtinv The inverse of the current time step.
 * \param approximate Use the approximate model showing bouncing.
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). The iterative method to solve the contact
 * problem is e.g. described in the article "A matrix-free cone complementarity approach for
 * solving large-scale, nonsmooth, rigid body dynamics" by A. Tasora and M. Anitescu in Computer
 * Methods in Applied Mechanics and Engineering (Volume 200, Issues 5–8, 15 January 2011,
 * Pages 439-453). The contact model is purely quadratic and convergence should be good but depends
 * on a parameter. The one-contact problem has a unique solution. The frictional reactions
 * for a dynamic contact converge to those that directly oppose slip. However, the contact is
 * not perfectly inelastic for dynamic contacts but bounces. These vertical motions tend to
 * go to zero for smaller time steps and can be interpreted as exaggerated vertical motions
 * coming from micro asperities (see "Optimization-based simulation of nonsmooth rigid multibody
 * dynamics" by M. Anitescu in Mathematical Programming (Volume 105, Issue 1, January 2006, Pages
 * 113-143). These motions can be prevented by a small change in the iteration proposed in "The
 * bipotential method: a constructive approach to design the complete contact law with friction and
 * improved numerical algorithms" by G. De Saxce and Z-Q. Feng in Mathematical and Computer
 * Modelling (Volume 28, Issue 4, 1998, Pages 225-245). Which iteration is used is controlled with
 * the approximate parameter.
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxInelasticCoulombContactsByOrthogonalProjections(size_t cid, real_t dtinv, bool approximate, CAccessor &ca, PAccessor& pa )
{
   real_t delta_max( 0 );

   // Remove velocity corrections of this contact's reaction.
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));

   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   Vec3 gdot    ( ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) -
                  ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) +
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) -
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) /* + diag_[i] * p */ );

   // Change from the global world frame to the contact frame
   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
   Vec3 gdot_nto( contactframe.getTranspose() * gdot );

   // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
   gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + ca.getDistance(cid) ) * dtinv;


   const real_t w( 1 ); // w > 0
   Vec3 p_cf( contactframe.getTranspose() * ca.getP(cid) );
   if( approximate ) {
      // Calculate next iterate (Anitescu/Tasora).
      p_cf = p_cf - w * ( ca.getDiag_nto(cid) * p_cf + gdot_nto );
   }
   else {
      // Calculate next iterate (De Saxce/Feng).
      Vec3 tmp( ca.getDiag_nto(cid) * p_cf + gdot_nto );
      tmp[0] += std::sqrt( math::sq( tmp[1] ) + math::sq( tmp[2] ) ) * ca.getMu(cid);
      p_cf = p_cf - w * tmp;
   }

   // Project.
   real_t flimit( ca.getMu(cid) * p_cf[0] );
   real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
   if( p_cf[0] > 0 && fsq < flimit * flimit ) {
      // Unconstrained minimum is in cone leading to a static contact and no projection
      // is necessary.
   }
   else if( p_cf[0] < 0 && fsq < math::sq( p_cf[0] )/ math::sq( ca.getMu(cid) ) ) {
      // Unconstrained minimum is in dual cone leading to a separating contact where no contact
      // reaction is present (the unconstrained minimum is projected to the tip of the cone).
      p_cf = Vec3();
   }
   else {
      // The contact is dynamic.
      real_t f( std::sqrt( fsq ) );
      p_cf[0] = ( f * ca.getMu(cid) + p_cf[0] ) / ( math::sq( ca.getMu(cid) ) + 1 );
      real_t factor( ca.getMu(cid) * p_cf[0] / f );
      p_cf[1] *= factor;
      p_cf[2] *= factor;
   }

   Vec3 p_wf( contactframe * p_cf );
   Vec3 dp( ca.getP(cid) - p_wf );
   delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

   ca.getPRef(cid) = p_wf;

   // Apply impulse right away.
   pa.getDvRef(bId1) += pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) += pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) -= pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) -= pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************

/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with the generalized maximum dissipation principle for friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Dynamic solutions are computed by
 * minimizing the kinetic energy along the intersection of the plane of maximum compression and
 * the friction cone.
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxInelasticGeneralizedMaximumDissipationContacts(size_t cid, real_t dtinv, CAccessor &ca, PAccessor& pa )
{
   real_t delta_max( 0 );

   // Remove velocity corrections of this contact's reaction.
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
   pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
   pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
   pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));

   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   Vec3 gdot    ( ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) -
                  ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) +
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) -
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) /* + diag_[i] * p */ );

   // Change from the global world frame to the contact frame
   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
   Vec3 gdot_nto( contactframe.getTranspose() * gdot );

   // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
   gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + ca.getDistance(cid) ) * dtinv;


   if( gdot_nto[0] >= 0 ) {
      // Contact is separating if no contact reaction is present at contact i.
      delta_max = std::max( delta_max, std::max( std::abs( ca.getP(cid)[0] ), std::max( std::abs( ca.getP(cid)[1] ), std::abs( ca.getP(cid)[2] ) ) ) );
      ca.getPRef(cid) = Vec3();

      // No need to apply zero impulse.
   }
   else {
      // Contact is persisting (either static or dynamic).

      // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
      Vec3 p_cf( -( ca.getDiag_nto_inv(cid) * gdot_nto ) );

      // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
      // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
      // b = [0.01 -1 -1]';
      // A\b    \approx [-0.19 -2.28 -1.21]'
      // eig(A) \approx [ 0.40  0.56  1.04]'

      real_t flimit( ca.getMu(cid) * p_cf[0] );
      real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
      if( fsq > flimit * flimit || p_cf[0] < 0 ) {
         // Contact cannot be static so it must be dynamic.
         // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

         // \breve{x}^0 = p_cf[1..2]

         // Eliminate normal component from 3x3 system: ca.getDiag_nto(cid)*p_cf + gdot_nto => \breve{A} \breve{x} - \breve{b}
         const real_t invA_nn( math::inv( ca.getDiag_nto(cid)(0, 0) ) );
         const real_t offdiag( ca.getDiag_nto(cid)(1, 2) - invA_nn * ca.getDiag_nto(cid)(0, 1) * ca.getDiag_nto(cid)(0, 2) );
         Mat2 A_breve( ca.getDiag_nto(cid)(1, 1) - invA_nn *math::sq( ca.getDiag_nto(cid)(0, 1) ), offdiag, offdiag, ca.getDiag_nto(cid)(2, 2) - invA_nn *math::sq( ca.getDiag_nto(cid)(0, 2) ) );
         Vec2 b_breve( -gdot_nto[1] + invA_nn * ca.getDiag_nto(cid)(0, 1) * gdot_nto[0], -gdot_nto[2] + invA_nn * ca.getDiag_nto(cid)(0, 2) * gdot_nto[0] );

         const real_t shiftI( std::atan2( -ca.getDiag_nto(cid)(0, 2), ca.getDiag_nto(cid)(0, 1) ) );
         const real_t shiftJ( std::atan2( -p_cf[2], p_cf[1] ) );
         const real_t a3( std::sqrt(math::sq( ca.getDiag_nto(cid)(0, 1) ) +math::sq( ca.getDiag_nto(cid)(0, 2) ) ) );
         const real_t fractionI( -ca.getDiag_nto(cid)(0, 0) / ( ca.getMu(cid) * a3 ) );
         const real_t fractionJ( std::min( invA_nn * ca.getMu(cid) * ( ( -gdot_nto[0] ) / std::sqrt( fsq ) - a3 * std::cos( shiftI - shiftJ ) ), real_c( 1 ) ) );

         // Search interval determination.
         real_t alpha_left, alpha_right;
         if( fractionJ < -1 ) {
            // J is complete
            const real_t angleI( std::acos( fractionI ) );
            alpha_left = -angleI - shiftI;
            alpha_right = +angleI - shiftI;
            if( alpha_left < 0 ) {
               alpha_left += 2 * math::pi;
               alpha_right += 2 * math::pi;
            }
         }
         else if( ca.getDiag_nto(cid)(0, 0) > ca.getMu(cid) * a3 ) {
            // I is complete
            const real_t angleJ( std::acos( fractionJ ) );
            alpha_left = -angleJ - shiftJ;
            alpha_right = +angleJ - shiftJ;
            if( alpha_left < 0 ) {
               alpha_left += 2 * math::pi;
               alpha_right += 2 * math::pi;
            }
         }
         else {
            // neither I nor J is complete
            const real_t angleJ( std::acos( fractionJ ) );
            real_t alpha1_left( -angleJ - shiftJ );
            real_t alpha1_right( +angleJ - shiftJ );
            if( alpha1_left < 0 ) {
               alpha1_left += 2 * math::pi;
               alpha1_right += 2 * math::pi;
            }
            const real_t angleI( std::acos( fractionI ) );
            real_t alpha2_left( -angleI - shiftI );
            real_t alpha2_right( +angleI - shiftI );
            if( alpha2_left < 0 ) {
               alpha2_left += 2 * math::pi;
               alpha2_right += 2 * math::pi;
            }

            // Swap intervals if second interval does not start right of the first interval.
            if( alpha1_left > alpha2_left ) {
               std::swap( alpha1_left, alpha2_left );
               std::swap( alpha1_right, alpha2_right );
            }

            if( alpha2_left > alpha1_right ) {
               alpha2_right -= 2*math::pi;
               if( alpha2_right > alpha1_right ) {
                  // [alpha1_left; alpha1_right] \subset [alpha2_left; alpha2_right]
               }
               else {
                  // [alpha2_left; alpha2_right] intersects the left end of [alpha1_left; alpha1_right]
                  alpha1_right = alpha2_right;
               }
            }
            else {
               alpha1_left = alpha2_left;
               if( alpha2_right > alpha1_right ) {
                  // [alpha2_left; alpha2_right] intersects the right end of [alpha1_left; alpha1_right]
               }
               else {
                  // [alpha2_left; alpha2_right] \subset [alpha1_left; alpha1_right]
                  alpha1_right = alpha2_right;
               }
            }

            alpha_left = alpha1_left;
            alpha_right = alpha1_right;
         }

         const real_t phi( real_c(0.5) * ( real_c(1) + std::sqrt( real_c( 5 ) ) ) );
         real_t alpha_mid( ( alpha_right + alpha_left * phi ) / ( 1 + phi ) );
         Vec2 x_mid;
         real_t f_mid;

         {
            real_t r_ub = ca.getMu(cid) * ( -gdot_nto[0] ) / ( ca.getDiag_nto(cid)(0, 0) + ca.getMu(cid) * a3 * std::cos( alpha_mid + shiftI ) );
            if( r_ub < 0 )
               r_ub = math::Limits<real_t>::inf();
            x_mid = Vec2( r_ub * std::cos( alpha_mid ), r_ub * std::sin( alpha_mid ) );
            f_mid = real_c(0.5) * x_mid * ( A_breve * x_mid ) - x_mid * b_breve;
         }

         bool leftlarger = false;
         for( size_t k = 0; k < maxSubIterations_; ++k ) {
            real_t alpha_next( alpha_left + ( alpha_right - alpha_mid ) );
            real_t r_ub = ca.getMu(cid) * ( -gdot_nto[0] ) / ( ca.getDiag_nto(cid)(0, 0) + ca.getMu(cid) * a3 * std::cos( alpha_next + shiftI ) );
            if( r_ub < 0 )
               r_ub = math::Limits<real_t>::inf();
            Vec2 x_next( r_ub * std::cos( alpha_next ), r_ub * std::sin( alpha_next ) );
            real_t f_next( real_c(0.5) * x_next * ( A_breve * x_next ) - x_next * b_breve );

            //WALBERLA_LOG_WARNING( "[(" << alpha_left << ", ?); (" << alpha_mid << ", " << f_mid << "); (" << alpha_right << ", ?)] <- (" << alpha_next << ", " << f_next << ")" );
            //WALBERLA_LOG_WARNING( "left: " << alpha_mid - alpha_left << "  right: " << alpha_right - alpha_mid << "  ll: " << leftlarger );
            //WALBERLA_ASSERT(leftlarger ? (alpha_mid - alpha_left > alpha_right - alpha_mid) : (alpha_mid - alpha_left < alpha_right - alpha_mid), "ll inconsistent!" );

            if (leftlarger) {
               // left interval larger
               if( f_next < f_mid ) {
                  alpha_right = alpha_mid;
                  alpha_mid   = alpha_next;
                  x_mid       = x_next;
                  f_mid       = f_next;
                  leftlarger = true;
               }
               else {
                  alpha_left  = alpha_next;
                  leftlarger = false;
               }
            }
            else {
               // right interval larger
               if( f_next < f_mid ) {
                  alpha_left = alpha_mid;
                  alpha_mid  = alpha_next;
                  x_mid      = x_next;
                  f_mid      = f_next;
                  leftlarger = false;
               }
               else {
                  alpha_right = alpha_next;
                  leftlarger = true;
               }
            }
         }
         //WALBERLA_LOG_DETAIL( "dalpha = " << alpha_right - alpha_left << "\n");
         {
            real_t alpha_init( std::atan2( p_cf[2], p_cf[1] ) );
            real_t r_ub = ca.getMu(cid) * ( -gdot_nto[0] ) / ( ca.getDiag_nto(cid)(0, 0) + ca.getMu(cid) * a3 * std::cos( alpha_init + shiftI ) );
            if( r_ub < 0 )
               r_ub = math::Limits<real_t>::inf();
            Vec2 x_init( r_ub * std::cos( alpha_init ), r_ub * std::sin( alpha_init ) );
            real_t f_init( real_c(0.5) * x_init * ( A_breve * x_init ) - x_init * b_breve );

            if( f_init < f_mid )
            {
               x_mid = x_init;
               WALBERLA_LOG_DETAIL( "Replacing solution by primitive dissipative solution (" << f_init << " < " << f_mid << " at " << alpha_init << " vs. " << alpha_mid << ").\n");
            }
         }

         p_cf[0] = invA_nn * ( -gdot_nto[0] - ca.getDiag_nto(cid)(0, 1) * x_mid[0] - ca.getDiag_nto(cid)(0, 2) * x_mid[1] );
         p_cf[1] = x_mid[0];
         p_cf[2] = x_mid[1];
         //WALBERLA_LOG_DETAIL( "Contact #" << i << " is dynamic." );
      }
      else {
         // Contact is static.
         //WALBERLA_LOG_DETAIL( "Contact #" << i << " is static." );
      }
      Vec3 p_wf( contactframe * p_cf );
      Vec3 dp( ca.getP(cid) - p_wf );
      delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );
      //WALBERLA_LOG_DETAIL( "Contact reaction in contact frame: " << p_cf << "\nContact action in contact frame: " << ca.getDiag_nto(cid)*p_cf + gdot_nto );

      ca.getPRef(cid) = p_wf;

      // Apply impulse right away.
      pa.getDvRef(bId1) += pa.getInvMass(bId1) * ca.getP(cid);
      pa.getDwRef(bId1) += pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
      pa.getDvRef(bId2) -= pa.getInvMass(bId2) * ca.getP(cid);
      pa.getDwRef(bId2) -= pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));
   }
   return delta_max;
}
//*************************************************************************************************

//*************************************************************************************************
//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is from the paper by Tschisgale et al.
 * "A constraint-based collision model for Cosserat rods" and works with a coefficient of
 * restitution cor with 0 < cor <= 1.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 */
template <typename CAccessor, typename PAccessor>
inline real_t HCSITSRelaxationStep::relaxInelasticContactsByProjectedGaussSeidel(size_t cid, real_t /*dtinv*/, CAccessor &ca, PAccessor& pa )
{
   real_t delta_max( 0 );
   size_t bId1 = ca.getId1(cid);
   size_t bId2 = ca.getId2(cid);
   // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
   // This is the negative contact velocity as for the other methods
   Vec3 gdot    ( ( pa.getLinearVelocity(bId2)     + pa.getDv(bId2) ) -
                  ( pa.getLinearVelocity(bId1)     + pa.getDv(bId1) ) +
                  ( pa.getAngularVelocity(bId2)    + pa.getDw(bId2) ) % ca.getR2(cid) -
                  ( pa.getAngularVelocity(bId1)    + pa.getDw(bId1) ) % ca.getR1(cid) /* + diag_[i] * p */ );

   // Use restitution hypothesis
   gdot = gdot + getCor()*(gdot*ca.getNormal(cid))*ca.getNormal(cid);

   Mat3 contactframe( ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );

   //
   // Turn system matrix back to the world frame.
   Mat3 K = contactframe * ca.getDiag_nto(cid) *  contactframe.getTranspose();
   Mat3 Kinv = contactframe * ca.getDiag_nto_inv(cid) * contactframe.getTranspose();

   // Compute impulse necessary in the world frame
   Vec3 p_wf( - Kinv * gdot );



   if( gdot * ca.getNormal(cid) <= 0 ) {
      // Contact is separating if no contact reaction is present at contact i.
      delta_max = std::max( delta_max, std::max( std::abs( ca.getP(cid)[0] ), std::max( std::abs( ca.getP(cid)[1] ), std::abs( ca.getP(cid)[2] ) ) ) );
      ca.getPRef(cid) = Vec3();
      // No need to apply zero impulse.
   }
   else {
      // Dissect the impuls in a tangential and normal directions
      // Use the inverted contact frame with -n as in the publication
      Mat3 reversedContactFrame( -ca.getNormal(cid), ca.getT(cid), ca.getO(cid) );
      Vec3 p_rcf( reversedContactFrame.getTranspose() * p_wf );

      // Check the frictional limit
      real_t flimit( ca.getMu(cid) * p_rcf[0] );

      // Square of tangential components of the impulse
      real_t fsq( p_rcf[1] * p_rcf[1] + p_rcf[2] * p_rcf[2] );

      if( fsq > flimit * flimit || p_rcf[0] < 0 ) {
         // Contact cannot be static so it must be dynamic.
         // Calculate corrected contact reaction by the method of Tschigale et al.
         // Normal component is close to 0 or negative, and we cannot asses a tangential direction
         // Treat this case with the 0 impulse
         if (fsq < real_t(1e-8)) {
            delta_max = std::max(delta_max, std::max(std::abs(ca.getP(cid)[0]),
                                                     std::max(std::abs(ca.getP(cid)[1]), std::abs(ca.getP(cid)[2]))));
            ca.getPRef(cid) = Vec3();
         } else {
            // tangential direction of sliding
            Vec3 tan_dir = ((p_rcf[1] * ca.getT(cid)) + (p_rcf[2] * ca.getO(cid))).getNormalized();
            // scalar magnitude of pv.
            real_t abspv = ((K * p_wf) * ca.getNormal(cid)) /
                           ((K * (-ca.getNormal(cid) + ca.getMu(cid) * tan_dir)) * ca.getNormal(cid));
            p_wf = abspv * (-ca.getNormal(cid) + ca.getMu(cid) * tan_dir);

         }
      }
      // Calculate variation
      Vec3 dp(ca.getP(cid) - p_wf);
      delta_max = std::max(delta_max, std::max(std::abs(dp[0]), std::max(std::abs(dp[1]), std::abs(dp[2]))));
      ca.getPRef(cid) = p_wf;

      // Apply impulse right away.
      pa.getDvRef(bId1) -= pa.getInvMass(bId1) * ca.getP(cid);
      pa.getDwRef(bId1) -= pa.getInvInertiaBF(bId1) * (ca.getR1(cid) % ca.getP(cid));
      pa.getDvRef(bId2) += pa.getInvMass(bId2) * ca.getP(cid);
      pa.getDwRef(bId2) += pa.getInvInertiaBF(bId2) * (ca.getR2(cid) % ca.getP(cid));
   }
   return delta_max;
}
//*************************************************************************************************


} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
