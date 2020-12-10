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
//! \file   DEMTangentialCollision.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/common/ParticleFunctions.h"

#include "mesa_pd/data/ParticleAccessor.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/SemiImplicitEuler.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include <iostream>

namespace walberla {
namespace mesa_pd {

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
 * Simulates oblique sphere-wall collision and checks rebound angle, i.e. the tangential part of the collision model.
 *
 */
int main( int argc, char ** argv )
{
   walberla::mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   real_t impactAngle = real_t(0);
   real_t dt         = real_c(2e-7);
   real_t frictionCoeff_s = real_t(0.8);
   real_t frictionCoeff_d = real_t(0.125);
   std::string filename = "TangentialCollision.txt";
   real_t collisionDuration = real_t(10);
   real_t nu = real_t(0.22); //Poissons ratio
   bool useVelocityVerlet = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--impactAngle" ) == 0 ) { impactAngle = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--dt" ) == 0 ) { dt = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--staticFriction" ) == 0 ) { frictionCoeff_s = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--dynamicFriction" ) == 0 ) { frictionCoeff_d = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--collisionDuration" ) == 0 ) { collisionDuration = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--nu" ) == 0 ) { nu = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--filename" ) == 0 ) { filename = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useVV" )   == 0 ) { useVelocityVerlet = true; continue; }
   }

   WALBERLA_LOG_INFO_ON_ROOT("******************************************************");
   WALBERLA_LOG_INFO_ON_ROOT("**                  NEW SIMULATION                  **");
   WALBERLA_LOG_INFO_ON_ROOT("******************************************************");
   WALBERLA_LOG_INFO_ON_ROOT("impactAngle = " << impactAngle);
   WALBERLA_LOG_INFO_ON_ROOT("dt = " << dt);
   WALBERLA_LOG_INFO_ON_ROOT("frictionCoeff_s = " << frictionCoeff_s);
   WALBERLA_LOG_INFO_ON_ROOT("frictionCoeff_d = " << frictionCoeff_d);
   WALBERLA_LOG_INFO_ON_ROOT("collisionDuration = " << collisionDuration);
   WALBERLA_LOG_INFO_ON_ROOT("nu = " << nu);

   real_t radius     = real_c(0.00159);
   real_t density    = real_c(2500);
   real_t restitutionCoeff = real_t(0.83);
   real_t collisionTime =  collisionDuration * dt;

   //init data structures
   auto ps = walberla::make_shared<data::ParticleStorage>(2);
   auto ss = walberla::make_shared<data::ShapeStorage>();
   using ParticleAccessor_T = ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   auto sphereShape = ss->create<data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia(density);

   const real_t particleMass =  real_t(1) / ss->shapes[sphereShape]->getInvMass();
   const real_t Mij = particleMass; // * particleMass / ( real_t(2) * particleMass ); // Mij = M for sphere-wall collision
   const real_t lnDryResCoeff = std::log(restitutionCoeff);

   // normal material aprameters
   const real_t stiffnessN = math::pi * math::pi * Mij / ( collisionTime * collisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ))  );
   const real_t dampingN = - real_t(2) * std::sqrt( Mij * stiffnessN ) *
   ( lnDryResCoeff / std::sqrt( math::pi * math::pi + ( lnDryResCoeff * lnDryResCoeff ) ) );

   WALBERLA_LOG_INFO_ON_ROOT("normal: stiffness = " << stiffnessN << ", damping = " << dampingN);

   const real_t lnDryResCoeffTangential = lnDryResCoeff; // std::log(0.31); //TODO: was same as in normal direction
   const real_t kappa = real_t(2) * ( real_t(1) - nu ) / ( real_t(2) - nu ) ;
   const real_t stiffnessT = kappa * Mij * math::pi * math::pi / ( collisionTime *  collisionTime );
   const real_t dampingT = real_t(2) * std::sqrt(Mij * stiffnessT) * ( - lnDryResCoeffTangential ) / ( std::sqrt( math::pi * math::pi + lnDryResCoeffTangential * lnDryResCoeffTangential ));

   WALBERLA_LOG_INFO_ON_ROOT("tangential: kappa = " << kappa << ", stiffness T = " << stiffnessT << ", damping T = " << dampingT);

   real_t uNin = real_t(1);
   real_t uTin = uNin * impactAngle;

   // create sphere
   data::Particle&& p = *ps->create();
   p.setPosition(Vec3(0,0,2*radius));
   p.setLinearVelocity(Vec3(uTin, 0., -uNin));
   p.setType(0);

   // create plane
   data::Particle&& p0 = *ps->create(true);
   p0.setPosition(Vec3(0,0,0));
   p0.setShapeID(ss->create<data::HalfSpace>(Vector3<real_t>(0,0,1)));
   p0.setType(0);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::FIXED);

   // velocity verlet
   kernel::VelocityVerletPreForceUpdate  vvPreForce( dt );
   kernel::VelocityVerletPostForceUpdate vvPostForce( dt );

   // explicit euler
   kernel::SemiImplicitEuler implEuler( dt );

   // collision response
   collision_detection::AnalyticContactDetection     acd;
   kernel::DoubleCast           double_cast;
   kernel::LinearSpringDashpot  dem(1);
   mpi::ReduceContactHistory    rch;
   dem.setStiffnessN(0,0,stiffnessN);
   dem.setStiffnessT(0,0,stiffnessT);
   dem.setDampingN(0,0,dampingN);
   dem.setDampingT(0,0,dampingT);
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
         //auto force = accessor->getForce(0);
         //WALBERLA_LOG_INFO(steps << ": penetration = " << penetration << " || vel = " << accessor->getLinearVelocity(0) << " || force = " << force);
      }
      rch(*ps);

      if(useVelocityVerlet) vvPostForce(0,*accessor);
      else implEuler(0, *accessor);

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
