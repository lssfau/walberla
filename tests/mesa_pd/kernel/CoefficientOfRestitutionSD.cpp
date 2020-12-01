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
//! \file CoefficientOfRestitutionSD.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"

#include "mesa_pd/data/ParticleAccessor.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/SemiImplicitEuler.h"
#include "mesa_pd/kernel/SpringDashpot.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include <string>

namespace dem_integrator_accuracy {

using namespace walberla;
using namespace walberla::mesa_pd;


class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss)
         : ParticleAccessor(ps)
         , ss_(ss)
   {}

   const auto& getInvMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)]->getInvMass();}

   const auto& getInvInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)]->getInvInertiaBF();}

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};


/*
 * Tests the integrator accuracy for a DEM simulation by comparing the given coefficient of restitution to the simulated one.
 * For that, the velocity after a single sphere-wall collision is divided by the initial velocity before the simulation.
 * The parameters of the DEM are chosen such as to (analytically) yield the desried coefficient of restitution.
 *
 * The simulation can be adapted via command line arguments.
 *
 * Currently compared integrators:
 *  - explicit euler (default)
 *  - velocity verlet (--useVV)
 */
int main( int argc, char** argv )
{
   mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   // parameters
   real_t radius = real_t(5);
   real_t dt = real_t(0.1);
   real_t restitutionCoeff = real_t(0.83);
   real_t densitySphere = real_t(1.5);
   real_t collisionTime = real_t(10);
   bool useVelocityVerlet = false;
   real_t uIn = real_t(0.1);

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--useVV" )   == 0 ) { useVelocityVerlet = true; continue; }
      if( std::strcmp( argv[i], "--dt" )      == 0 ) { dt = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Tc" )      == 0 ) { collisionTime = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--radius" )  == 0 ) { radius = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--e" )       == 0 ) { restitutionCoeff = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--density" ) == 0 ) { densitySphere = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--uIn" )     == 0 ) { uIn = real_c(std::atof( argv[++i] )); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   //init data structures
   auto ps = walberla::make_shared<data::ParticleStorage>(2);
   auto ss = walberla::make_shared<data::ShapeStorage>();
   using ParticleAccessor_T = ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   auto sphereShape = ss->create<data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia(densitySphere);

   const real_t particleMass = real_t(1) / ss->shapes[sphereShape]->getInvMass();
   const real_t Mij = particleMass; // Mij = M for sphere-wall collision
   const real_t lnDryResCoeff = std::log(restitutionCoeff);
   const real_t stiffnessN = math::pi * math::pi * Mij / ( collisionTime * collisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ))  );
   const real_t dampingN = - real_t(2) * std::sqrt( Mij * stiffnessN ) * ( lnDryResCoeff / std::sqrt( math::pi * math::pi + ( lnDryResCoeff * lnDryResCoeff ) ) );

   WALBERLA_LOG_INFO("dt = " << dt << ", Tc = " << collisionTime << ", coefficient of restitution = " << restitutionCoeff);
   WALBERLA_LOG_INFO(" -> mass " << particleMass << ", collision duration = " << collisionTime / dt << ", stiffness = " << stiffnessN << ", damping = " << dampingN);

   // create sphere
   auto pos = Vec3(0,0,radius);
   auto linVel = Vec3(0,0,-uIn);
   data::Particle&& p = *ps->create();
   p.setPosition(pos);
   p.setLinearVelocity(linVel);
   p.setShapeID(sphereShape);
   p.setType(0);
   p.setForce(Vec3(0.));
   p.setOldForce(Vec3(0.));

   // create plane
   data::Particle&& p0 = *ps->create(true);
   p0.setPosition(Vec3(0,0,0));
   p0.setShapeID(ss->create<data::HalfSpace>( Vector3<real_t>(0,0,1)) );
   p0.setType(0);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::FIXED);

   collision_detection::AnalyticContactDetection acd;
   kernel::DoubleCast       double_cast;

   kernel::SemiImplicitEuler implEuler(dt);
   kernel::VelocityVerletPreForceUpdate  vvPreForce( dt );
   kernel::VelocityVerletPostForceUpdate vvPostForce( dt );

   kernel::SpringDashpot dem(1);
   dem.setStiffness(0,0,stiffnessN);
   dem.setDampingN(0,0,dampingN);

   uint_t steps = 0;
   real_t maxPenetration = real_t(0);
   do
   {
      if(useVelocityVerlet) vvPreForce(0,*accessor);

      if (double_cast(0, 1, *accessor, acd, *accessor ))
      {
         real_t penetration = acd.getPenetrationDepth();
         maxPenetration = std::max( maxPenetration, std::abs(penetration));

         dem(acd.getIdx1(), acd.getIdx2(), *accessor, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
      }
      auto force = accessor->getForce(0);

      if(useVelocityVerlet) vvPostForce(0,*accessor);
      else implEuler(0, *accessor);

      WALBERLA_LOG_INFO(steps << ": penetration = " << acd.getPenetrationDepth() << " || vel = " << accessor->getLinearVelocity(0)[2] << " || force = " << force[2]);

      ++steps;
   } while (double_cast(0, 1, *accessor, acd, *accessor ));

   real_t simulatedCoefficientOfRestitution = -accessor->getLinearVelocity(0)[2] / linVel[2];
   real_t relativeError = ( simulatedCoefficientOfRestitution - restitutionCoeff ) / restitutionCoeff;
   WALBERLA_LOG_INFO("coefficient of restitution = " << simulatedCoefficientOfRestitution << " -> error = " << relativeError * 100. << "%");
   WALBERLA_LOG_INFO("collision steps = " << steps);
   WALBERLA_LOG_INFO("Max penetration = " << maxPenetration << " -> " << maxPenetration / radius * 100. << "% of radius");

   if( useVelocityVerlet )
   {
      WALBERLA_CHECK_LESS(relativeError, real_t(0.01), "Error in simulated coefficient of restitution too large: " << simulatedCoefficientOfRestitution << " vs ref " << restitutionCoeff);
   }
   else
   {
      WALBERLA_CHECK_LESS(relativeError, real_t(0.03), "Error in simulated coefficient of restitution too large: " << simulatedCoefficientOfRestitution << " vs ref " << restitutionCoeff);
   }


   return EXIT_SUCCESS;
}
} // namespace dem_integrator_accuracy

int main( int argc, char** argv )
{
   return dem_integrator_accuracy::main(argc, argv);
}
