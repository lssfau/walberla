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
//! \file   HydrodynamicForceOnMultipleBlocks.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"


namespace hydrodynamic_force_blocks
{

using namespace walberla;


template< typename ParticleAccessor_T>
class SpherePropertyEvaluator
{
public:
   SpherePropertyEvaluator( const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid) :
         ac_( ac ), sphereUid_( sphereUid )
   {  }

   void operator()()
   {
      Vector3<real_t> pos(real_t(0));
      Vector3<real_t> transVel(real_t(0));
      Vector3<real_t> angularVel(real_t(0));
      Vector3<real_t> force(real_t(0));
      Vector3<real_t> torque(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos = ac_->getPosition(idx);
            transVel = ac_->getLinearVelocity(idx);
            angularVel = ac_->getAngularVelocity(idx);
            force = ac_->getForce(idx);
            torque = ac_->getTorque(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos[0], mpi::SUM );
         mpi::allReduceInplace( pos[1], mpi::SUM );
         mpi::allReduceInplace( pos[2], mpi::SUM );

         mpi::allReduceInplace( transVel[0], mpi::SUM );
         mpi::allReduceInplace( transVel[1], mpi::SUM );
         mpi::allReduceInplace( transVel[2], mpi::SUM );

         mpi::allReduceInplace( angularVel[0], mpi::SUM );
         mpi::allReduceInplace( angularVel[1], mpi::SUM );
         mpi::allReduceInplace( angularVel[2], mpi::SUM );

         mpi::allReduceInplace( force[0], mpi::SUM );
         mpi::allReduceInplace( force[1], mpi::SUM );
         mpi::allReduceInplace( force[2], mpi::SUM );

         mpi::allReduceInplace( torque[0], mpi::SUM );
         mpi::allReduceInplace( torque[1], mpi::SUM );
         mpi::allReduceInplace( torque[2], mpi::SUM );
      }

      position_ = pos;
      linearVel_ = transVel;
      angularVel_ = angularVel;
      force_ = force;
      torque_ = torque;
   }

   Vector3<real_t> getPosition() const { return position_; }
   Vector3<real_t> getLinearVel() const { return linearVel_; }
   Vector3<real_t> getAngularVel() const { return angularVel_; }
   Vector3<real_t> getForce() const { return force_; }
   Vector3<real_t> getTorque() const { return torque_; }


private:

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;

   Vector3<real_t> position_;
   Vector3<real_t> linearVel_;
   Vector3<real_t> angularVel_;
   Vector3<real_t> force_;
   Vector3<real_t> torque_;
};

template< typename ParticleAccessor_T>
void applyHydrodynamicForceTorqueOnSphere(ParticleAccessor_T & accessor, walberla::id_t sphereUid,
                                          Vector3<real_t> hydForce, Vector3<real_t> hydTorque)
{

   uint_t numberOfProcessesWithKnowledgeOfThisSphere = uint_t(0);
   size_t idx = accessor.uidToIdx(sphereUid);
   if( idx != accessor.getInvalidIdx())
   {
      ++numberOfProcessesWithKnowledgeOfThisSphere;
   }

   WALBERLA_MPI_SECTION() {
      mpi::allReduceInplace(numberOfProcessesWithKnowledgeOfThisSphere, mpi::SUM);
   }

   if( idx != accessor.getInvalidIdx())
   {
      accessor.setHydrodynamicForce(idx, hydForce / real_t(numberOfProcessesWithKnowledgeOfThisSphere) );
      accessor.setHydrodynamicTorque(idx, hydTorque / real_t(numberOfProcessesWithKnowledgeOfThisSphere) );
   }
}


/*
 * Two spheres travelling through several blocks. moved by constant (distributed) hydrodynamic force.
 * This checks setting hydrodynamic forces onto particles and force averaging in parallel setup.
 * Velocity of both spheres have to be equal throughout the simulation.
 *
 */
int main( int argc, char ** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   const Vector3<uint_t> domainSize( 4 * 32, 32, 32 );
   const Vector3<uint_t> numberOfBlocksPerDirection(4, 1, 1 );
   Vector3<uint_t> cellsPerBlockPerDirection( domainSize[0] / numberOfBlocksPerDirection[0],
                                              domainSize[1] / numberOfBlocksPerDirection[1],
                                              domainSize[2] / numberOfBlocksPerDirection[2] );

   real_t xPos1 = real_t(0);
   real_t xPos2 = real_t(45.36281); // random

   real_t radius = real_t(2);

   Vector3<real_t> hydForce( real_t(0.75), 0, 0);
   Vector3<real_t> hydTorque( real_t(0.75), 0, 0);

   real_t dt = real_t(1);
   real_t dx = real_t(1);
   uint_t timesteps = 400;

   bool averageForceTorqueOverTwoTimeSteps = true;
   bool useVelocityVerlet = false;


   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--timesteps" )           == 0 ) { timesteps = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" )    == 0 ) { averageForceTorqueOverTwoTimeSteps = false; continue; }
      if( std::strcmp( argv[i], "--useVV" )               == 0 ) { useVelocityVerlet = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   // setup as in a coupled simulation even though not needed here for mesa_pd only setup
   auto blocks = blockforest::createUniformBlockGrid( numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1], numberOfBlocksPerDirection[2],
                                                      cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx,
                                                      0, false, false,
                                                      true, false, false, //periodicity
                                                      false );


   auto rpdDomain = walberla::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(2);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia(real_t(1));

   Vector3<real_t> initialPosition1( xPos1, real_t(0.5) * real_c(domainSize[1]), real_t(0.5) * real_c(domainSize[2]));
   Vector3<real_t> initialPosition2( xPos2, real_t(0.5) * real_c(domainSize[1]), real_t(0.5) * real_c(domainSize[2]));

   walberla::id_t sphereUid1 = 0;
   if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), initialPosition1 ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(initialPosition1);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid1 = p.getUid();
   }
   mpi::allReduceInplace(sphereUid1, mpi::SUM);

   walberla::id_t sphereUid2 = 0;
   if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), initialPosition2 ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(initialPosition2);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid2 = p.getUid();
   }
   mpi::allReduceInplace(sphereUid2, mpi::SUM);


   // mesa_pd functionality
   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce( dt );
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce( dt );
   mesa_pd::kernel::ExplicitEuler explEulerIntegrator( dt );
   mesa_pd::mpi::ReduceProperty reduceProperty;

   std::function<void(void)> syncCall = [ps,rpdDomain](){
      const real_t overlap = real_t( 1.5 );
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   // coupling functionality
   lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;

   // evaluation functionality
   SpherePropertyEvaluator<ParticleAccessor_T> sphere1Eval(accessor, sphereUid1);
   SpherePropertyEvaluator<ParticleAccessor_T> sphere2Eval(accessor, sphereUid2);

   bool useOpenMP = false;
   for(uint_t t = 0; t < timesteps; ++t) {

      // set hydrodynamic force/torque distributed
      applyHydrodynamicForceTorqueOnSphere(*accessor, sphereUid1, hydForce, hydTorque);
      applyHydrodynamicForceTorqueOnSphere(*accessor, sphereUid2, hydForce, hydTorque);

      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);

      if (averageForceTorqueOverTwoTimeSteps) {
         if (t == 0) {
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor);
         }
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor);
      }

      if (useVelocityVerlet) {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         syncCall();
      }

      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor);

      reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

      // check sphere properties
      sphere1Eval();
      sphere2Eval();

      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getLinearVel()[0], sphere2Eval.getLinearVel()[0], "time step " << t << " mismatch in linear vel[0]");
      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getLinearVel()[1], sphere2Eval.getLinearVel()[1], "time step " << t << " mismatch in linear vel[1]");
      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getLinearVel()[2], sphere2Eval.getLinearVel()[2], "time step " << t << " mismatch in linear vel[2]");

      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getAngularVel()[0], sphere2Eval.getAngularVel()[0], "time step " << t << " mismatch in angular vel[0]");
      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getAngularVel()[1], sphere2Eval.getAngularVel()[1], "time step " << t << " mismatch in angular vel[1]");
      WALBERLA_CHECK_FLOAT_EQUAL(sphere1Eval.getAngularVel()[2], sphere2Eval.getAngularVel()[2], "time step " << t << " mismatch in angular vel[2]");

      // particle integration
      if( useVelocityVerlet ) ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
      else ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explEulerIntegrator, *accessor);

      syncCall();

      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

      //WALBERLA_LOG_INFO_ON_ROOT(t << " " << sphere1Eval.getPosition()[0] << " " << sphere1Eval.getLinearVel()[0] << " " << sphere1Eval.getForce()[0] << " " << sphere1Eval.getAngularVel()[0] << " " << sphere1Eval.getTorque()[0] );
   }

   return EXIT_SUCCESS;
}

} //namespace hydrodynamic_force_blocks

int main( int argc, char ** argv )
{
   return hydrodynamic_force_blocks::main(argc, argv);
}
