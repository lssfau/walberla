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
//! \file   ObliqueDryCollision.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/common/ParticleFunctions.h"

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include <iostream>

namespace walberla {
namespace mesa_pd {

/*
 * Simulates oblique sphere-wall collision and checks rebound angle, i.e. the tangential part of the collision model.
 *
 */
int main( int argc, char ** argv )
{
   walberla::mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   real_t uNin_SI = real_t(1.7); // m/s
   real_t diameter_SI = real_t(0.00318); // m
   //real_t density_SI = real_t(2500); // kg/m**3, not used

   real_t uNin       = real_t(0.02);
   real_t diameter   = real_t(20);
   real_t radius     = real_t(0.5) * diameter;
   real_t density    = real_c(2.5);

   // these values have actually no influence here and are just computed for completeness
   real_t dx_SI = diameter_SI / diameter;
   real_t dt_SI = uNin / uNin_SI * dx_SI;

   real_t impactAngle       = real_t(0);
   real_t dt                = real_t(0.1); // = (1 / #sub steps)
   real_t frictionCoeff_s   = real_t(0.8); // no influence
   real_t frictionCoeff_d   = real_t(0.12); // paper: 0.125+-0.007
   std::string filename     = "TangentialCollision.txt";
   real_t collisionTime     = real_t(80);
   real_t nu                = real_t(0.22); //Poissons ratio
   bool useVelocityVerlet   = false;
   real_t restitutionCoeff  = real_t(0.83);

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--impactAngle" ) == 0 ) { impactAngle = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--dt" ) == 0 ) { dt = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--staticFriction" ) == 0 ) { frictionCoeff_s = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--dynamicFriction" ) == 0 ) { frictionCoeff_d = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--collisionTime" ) == 0 ) { collisionTime = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--nu" ) == 0 ) { nu = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--filename" ) == 0 ) { filename = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useVV" )   == 0 ) { useVelocityVerlet = true; continue; }
      if( std::strcmp( argv[i], "--coefficientOfRestitution" ) == 0 ) { restitutionCoeff = real_c( std::atof( argv[++i] ) ); continue; }
   }

   WALBERLA_LOG_INFO_ON_ROOT("******************************************************");
   WALBERLA_LOG_INFO_ON_ROOT("**                  NEW SIMULATION                  **");
   WALBERLA_LOG_INFO_ON_ROOT("******************************************************");
   WALBERLA_LOG_INFO_ON_ROOT("dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT("dx_SI = " << dx_SI);
   WALBERLA_LOG_INFO_ON_ROOT("impactAngle = " << impactAngle);
   WALBERLA_LOG_INFO_ON_ROOT("dt = " << dt);
   WALBERLA_LOG_INFO_ON_ROOT("frictionCoeff_s = " << frictionCoeff_s);
   WALBERLA_LOG_INFO_ON_ROOT("frictionCoeff_d = " << frictionCoeff_d);
   WALBERLA_LOG_INFO_ON_ROOT("collisionTime = " << collisionTime);
   WALBERLA_LOG_INFO_ON_ROOT("nu = " << nu);
   WALBERLA_LOG_INFO_ON_ROOT("integrator = " << (useVelocityVerlet ? "Velocity Verlet" : "Explicit Euler"));
   WALBERLA_LOG_INFO_ON_ROOT("restitutionCoeff = " << restitutionCoeff);


   //init data structures
   auto ps = walberla::make_shared<data::ParticleStorage>(2);
   auto ss = walberla::make_shared<data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   auto sphereShape = ss->create<data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia(density);

   const real_t particleMass =  real_t(1) / ss->shapes[sphereShape]->getInvMass();
   const real_t Mij = particleMass; // * particleMass / ( real_t(2) * particleMass ); // Mij = M for sphere-wall collision
   const real_t kappa = real_t(2) * ( real_t(1) - nu ) / ( real_t(2) - nu ) ;    // from Thornton et al

   real_t uTin = uNin * impactAngle;

   // create sphere
   data::Particle&& p = *ps->create();
   p.setPosition(Vec3(0,0,2*radius));
   p.setLinearVelocity(Vec3(uTin, 0., -uNin));
   p.setInteractionRadius(radius);
   p.setType(0);

   // create plane
   data::Particle&& p0 = *ps->create(true);
   p0.setPosition(Vec3(0,0,0));
   p0.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p0.setShapeID(ss->create<data::HalfSpace>(Vector3<real_t>(0,0,1)));
   p0.setType(0);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::FIXED);

   // velocity verlet
   kernel::VelocityVerletPreForceUpdate  vvPreForce( dt );
   kernel::VelocityVerletPostForceUpdate vvPostForce( dt );

   // explicit euler
   kernel::ExplicitEuler explEuler( dt );

   // collision response
   collision_detection::AnalyticContactDetection acd;
   kernel::DoubleCast           double_cast;
   kernel::LinearSpringDashpot  dem(1);
   mpi::ReduceContactHistory    rch;
   dem.setStiffnessAndDamping(0,0,restitutionCoeff,collisionTime,kappa,Mij);
   dem.setFrictionCoefficientStatic(0,0,frictionCoeff_s);
   dem.setFrictionCoefficientDynamic(0,0,frictionCoeff_d);

   WALBERLA_LOG_DEVEL("begin: vel = " << p.getLinearVelocity() << ", contact vel: " << getVelocityAtWFPoint(0,*accessor,p.getPosition() + Vec3(0,0,-radius)) );

   uint_t steps = 0;
   real_t maxPenetration = real_t(0);
   do
   {
      if(useVelocityVerlet) vvPreForce(0,*accessor);

      real_t penetration;

      if (double_cast(0, 1, *accessor, acd, *accessor ))
      {
         penetration = acd.getPenetrationDepth();
         maxPenetration = std::max( maxPenetration, std::abs(penetration));

         dem(acd.getIdx1(), acd.getIdx2(), *accessor, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth(), dt);
         auto force = accessor->getForce(0);
         auto torque = accessor->getTorque(0);
         WALBERLA_LOG_INFO(steps << ": penetration = " << penetration << " || vel = " << accessor->getLinearVelocity(0) << " || force = " << force << ", torque = " << torque);
      }
      rch(*ps);

      if(useVelocityVerlet) vvPostForce(0,*accessor);
      else explEuler(0, *accessor);

      ++steps;
   } while (double_cast(0, 1, *accessor, acd, *accessor ) || p.getLinearVelocity()[2] < 0);

   real_t uTout = p.getLinearVelocity()[0] - radius * p.getAngularVelocity()[1];
   WALBERLA_LOG_DEVEL("end: linear vel = " << p.getLinearVelocity() << ", angular vel = " << p.getAngularVelocity());

   real_t reboundAngle = uTout / uNin;
   WALBERLA_LOG_INFO_ON_ROOT("gamma_in = " << impactAngle);
   WALBERLA_LOG_INFO_ON_ROOT("gamma_out = " << reboundAngle);

   WALBERLA_LOG_INFO_ON_ROOT("Thornton: sliding should occur if " << real_t(2) * impactAngle / ( frictionCoeff_d * ( real_t(1) + restitutionCoeff)) << " >= " << real_t(7) - real_t(1) / kappa );
   WALBERLA_LOG_INFO_ON_ROOT("Max penetration = " << maxPenetration << " -> " << maxPenetration / radius * 100. << "% of radius");

   std::ofstream file;
   file.open( filename.c_str(), std::ios::out | std::ios::app );
   file << impactAngle << " " << reboundAngle << "\n";
   file.close();

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
