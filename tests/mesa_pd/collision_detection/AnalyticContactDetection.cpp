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
//! \file   AnalyticContactDetection.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/kernel/DoubleCast.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <cmath>
#include <iostream>
#include <memory>

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

void checkRotation( const Vec3& from, const Vec3& to )
{
   WALBERLA_LOG_DEVEL_VAR(from);
   WALBERLA_LOG_DEVEL_VAR(to);
   real_t cos = dot(from, to) / from.length() / to.length();
   real_t angle = std::acos(cos);
   WALBERLA_LOG_DEVEL_VAR(angle);
   Vec3 axis = cross(from, to);
   Rot3 rot;
   rot.rotate(axis, angle);
   WALBERLA_LOG_DEVEL_VAR( rot );

   auto res = rot.getMatrix() * from;

   WALBERLA_LOG_DEVEL_VAR(res);

   WALBERLA_CHECK_FLOAT_EQUAL( res, to );
}

void checkSphereSphereCollision( )
{
   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape ac(ps, ss);

   //initialize particles
   const real_t radius  = real_t(0.5);
   auto smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   data::Particle&& sp0        = *ps->create();
   sp0.setPosition( Vec3(real_t(0), real_t(0), real_t(0)) );
   sp0.setShapeID(  smallSphere );

   auto dir = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   real_t shift = real_c(0.75);
   data::Particle&& sp1        = *ps->create();
   sp1.setPosition( dir * shift );
   sp1.setShapeID(  smallSphere );

   collision_detection::AnalyticContactDetection         acd;
   kernel::DoubleCast               double_cast;

   bool isInContact = false;
   isInContact = double_cast(0, 1, ac, acd, ac );

   //check two spheres in contact
   WALBERLA_CHECK( isInContact );
   WALBERLA_CHECK_EQUAL( acd.getIdx1(), 0);
   WALBERLA_CHECK_EQUAL( acd.getIdx2(), 1);
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactPoint(), dir * shift * real_t(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactNormal(), -dir );
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getPenetrationDepth(), shift - real_t(1) );

   //check order invariance
   collision_detection::AnalyticContactDetection         acd2;
   WALBERLA_CHECK( double_cast(1, 0, ac, acd2, ac ) );
   WALBERLA_CHECK_EQUAL( acd.getIdx1(), acd2.getIdx1());
   WALBERLA_CHECK_EQUAL( acd.getIdx2(), acd2.getIdx2());
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactPoint(), acd2.getContactPoint() );
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactNormal(), acd2.getContactNormal() );
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getPenetrationDepth(), acd2.getPenetrationDepth() );

   //no collision
   sp1.setPosition( dir * real_t(1.1) );
   WALBERLA_CHECK( !double_cast(1, 0, ac, acd, ac ) );
}

void checkSphereHalfSpaceCollision( )
{
   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape ac(ps, ss);

   //initialize particles
   const real_t radius  = real_t(1.0);
   auto smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   data::Particle&& sp0        = *ps->create();
   sp0.setPosition( Vec3(real_t(0), real_t(0), real_t(0)) );
   sp0.setShapeID(  smallSphere );

   auto dir = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   real_t shift = real_c(0.75);
   auto p1              = ps->create(true);
   p1->setPosition( dir * shift );
   p1->setShapeID(  ss->create<data::HalfSpace>( -dir ) );

   collision_detection::AnalyticContactDetection         acd;
   kernel::DoubleCast               double_cast;

   bool isInContact = false;
   isInContact = double_cast(0, 1, ac, acd, ac );

   //check sphere - half space contact
   WALBERLA_CHECK( isInContact );
   WALBERLA_CHECK_EQUAL( acd.getIdx1(), 0);
   WALBERLA_CHECK_EQUAL( acd.getIdx2(), 1);
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactPoint(), dir * shift * real_t(1.0));
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactNormal(), -dir );
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getPenetrationDepth(), shift - real_t(1) );

   auto pos = Vec3(shift, real_t(0), real_t(0));
   p1->setPosition(pos);
   p1->setShapeID(  ss->create<data::HalfSpace>( -Vec3(real_t(1), real_t(0), real_t(0)) ) );
   isInContact = double_cast(0, 1, ac, acd, ac );
   WALBERLA_CHECK( isInContact );
   WALBERLA_CHECK_EQUAL( acd.getIdx1(), 0);
   WALBERLA_CHECK_EQUAL( acd.getIdx2(), 1);
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactPoint(), pos);
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getContactNormal(), -Vec3(real_t(1), real_t(0), real_t(0)) );
   WALBERLA_CHECK_FLOAT_EQUAL( acd.getPenetrationDepth(), pos[0] - real_t(1) );
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla;
   using namespace walberla::mesa_pd;
   walberla::mesa_pd::checkRotation( Vec3(real_t(1), real_t(0), real_t(0)),
                                     Vec3(real_t(0), real_t(1), real_t(0)) );
   walberla::mesa_pd::checkSphereSphereCollision();
   walberla::mesa_pd::checkSphereHalfSpaceCollision();

   return EXIT_SUCCESS;
}
