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
//! \file   SpringDashpot.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/SpringDashpot.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

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

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();

   auto smallSphere = ss->create<data::Sphere>( real_t(2) );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2500));

   ParticleAccessorWithShape ac(ps, ss);

   data::Particle&& p1 = *ps->create();
   p1.getPositionRef() = Vec3(0,0,0);
   //p1.getAngularVelocityRef() = Vec3(1,1,-1).getNormalized();
   p1.getShapeIDRef()  = smallSphere;
   p1.getTypeRef()     = 0;

   data::Particle&& p2 = *ps->create();
   p2.getPositionRef() = Vec3(2,2,2);
   p2.getShapeIDRef()  = smallSphere;
   p2.getTypeRef()     = 0;

   // Init kernels
   kernel::SpringDashpot sd(1);
   sd.setStiffness(0, 0, real_t(8.11e6));
   sd.setDampingN (0, 0, real_t(6.86e1));
   sd.setDampingT (0, 0, real_t(6.86e1));
   sd.setFriction (0, 0, real_t(1.2));

   collision_detection::AnalyticContactDetection contact;
   kernel::DoubleCast       double_cast;
   WALBERLA_CHECK(double_cast(0, 1, ac, contact, ac ));

   // single contact test
   sd(contact.getIdx1(),
      contact.getIdx2(),
      ac,
      contact.getContactPoint(),
      contact.getContactNormal(),
      contact.getPenetrationDepth());
   std::for_each(ps->begin(), ps->end(), [](data::Particle&& p){ WALBERLA_LOG_DEVEL(p); });
   WALBERLA_CHECK_FLOAT_EQUAL( ps->getForce(0), -ps->getForce(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( ps->getForce(0), Vec3(1,1,1).getNormalized() * ((std::sqrt(real_t(12)) - 4) * sd.getStiffness(0, 0)) );

   // thread safety test
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for (int i = 0; i < 100; ++i)
      sd(contact.getIdx1(),
         contact.getIdx2(),
         ac,
         contact.getContactPoint(),
         contact.getContactNormal(),
         contact.getPenetrationDepth());

   WALBERLA_CHECK_FLOAT_EQUAL( ps->getForce(0), -ps->getForce(1) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( ps->getForce(0),
                                       real_t(101) * Vec3(1,1,1).getNormalized() * ((std::sqrt(real_t(12)) - 4) * sd.getStiffness(0, 0)),
                                       real_t(1e-6) );

   auto cor  = real_t(0.87);
   auto ct   = real_t(0.17);
   auto meff = real_t(0.65);
   sd.setParametersFromCOR(0, 0, cor, ct, meff);
   //WALBERLA_CHECK_FLOAT_EQUAL(sd.getStiffness(0,0), (math::pi*math::pi - std::log(cor)*std::log(cor)) / (ct*ct) * meff);
   //WALBERLA_CHECK_FLOAT_EQUAL(sd.getDampingN(0,0),  -real_t(2)*std::log(cor)/ct*meff);
   WALBERLA_CHECK_FLOAT_EQUAL(sd.calcCoefficientOfRestitution(0, 0, meff), cor);
   WALBERLA_CHECK_FLOAT_EQUAL(sd.calcCollisionTime(0, 0, meff), ct);

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
