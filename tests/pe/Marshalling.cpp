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
//! \file Marshalling.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Squirmer.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/Ellipsoid.h"
#include "pe/communication/rigidbody/Squirmer.h"
#include "pe/communication/DynamicMarshalling.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Materials.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::pe::communication;

using UnionT = Union<Sphere>;
using UnionID = UnionT *;
using UnionPtr = std::unique_ptr<UnionT>;

typedef std::tuple<Box, Capsule, Sphere, Squirmer, UnionT, Ellipsoid> BodyTuple ;

void testBox()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testBox ***");

   MaterialID iron = Material::find("iron");

   Box b1(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), Vec3(1,2,3), iron, false, true, false);
   b1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   b1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, b1);
   mpi::RecvBuffer rb(sb);

   auto bPtr = UnmarshalDynamically<BodyTuple>::execute(rb, Box::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   BoxID b2 = static_cast<BoxID>(bPtr.get());

   WALBERLA_CHECK_FLOAT_EQUAL(b1.getPosition(), b2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getLinearVel(), b2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getAngularVel(), b2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getLengths(), b2->getLengths());
   WALBERLA_CHECK_EQUAL(b1.getID(), b2->getID());
   WALBERLA_CHECK_EQUAL(b1.getSystemID(), b2->getSystemID());
}

void testCapsule()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testCapsule ***");

   MaterialID iron = Material::find("iron");

   Capsule c1(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), 5, 7, iron, false, false, false);
   c1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   c1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, c1);
   mpi::RecvBuffer rb(sb);

   auto cPtr = UnmarshalDynamically<BodyTuple>::execute(rb, Capsule::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   CapsuleID c2 = static_cast<CapsuleID> (cPtr.get());

   WALBERLA_CHECK_FLOAT_EQUAL(c1.getPosition(), c2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(c1.getLinearVel(), c2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(c1.getAngularVel(), c2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(c1.getRadius(), c2->getRadius());
   WALBERLA_CHECK_FLOAT_EQUAL(c1.getLength(), c2->getLength());
   WALBERLA_CHECK_EQUAL(c1.getID(), c2->getID());
   WALBERLA_CHECK_EQUAL(c1.getSystemID(), c2->getSystemID());
}

void testSphere()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testSphere ***");

   MaterialID iron = Material::find("iron");

   Sphere s1(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)),  Quat(), 5, iron, false, false, false);
   s1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   s1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, s1);
   mpi::RecvBuffer rb(sb);

   auto sPtr = UnmarshalDynamically<BodyTuple>::execute(rb, Sphere::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   SphereID s2 = static_cast<SphereID> (sPtr.get());

   WALBERLA_CHECK_FLOAT_EQUAL(s1.getPosition(), s2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(s1.getLinearVel(), s2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(s1.getAngularVel(), s2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(s1.getRadius(), s2->getRadius());
   WALBERLA_CHECK_EQUAL(s1.getID(), s2->getID());
   WALBERLA_CHECK_EQUAL(s1.getSystemID(), s2->getSystemID());
}

void testSquirmer()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testSquirmer ***");

   MaterialID iron = Material::find("iron");

   Squirmer s1(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), real_c(5), real_c(0.1), real_c(4.93), iron, false, false, false);

   mpi::SendBuffer sb;
   WALBERLA_ASSERT_UNEQUAL(Sphere::getStaticTypeID(), Squirmer::getStaticTypeID(), "Squirmer did not get its own type ID");
   WALBERLA_ASSERT_UNEQUAL(Sphere::getStaticTypeID(), s1.getStaticTypeID(), "Squirmer did not get its own type ID");
   MarshalDynamically<BodyTuple>::execute(sb, s1);
   mpi::RecvBuffer rb(sb);

   auto sPtr = UnmarshalDynamically<BodyTuple>::execute(rb, Squirmer::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   SquirmerID s2 = static_cast<SquirmerID> (sPtr.get());

   WALBERLA_CHECK_FLOAT_EQUAL(s1.getSquirmerVelocity(), s2->getSquirmerVelocity());
   WALBERLA_CHECK_FLOAT_EQUAL(s1.getSquirmerBeta(), s2->getSquirmerBeta());
}

void testEllipsoid()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testEllipsoid ***");

   MaterialID iron = Material::find("iron");

   Ellipsoid e1(759847, 1234795, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), Vec3(3,1,5), iron, false, false, false);
   e1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   e1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, e1);
   mpi::RecvBuffer rb(sb);

   auto ePtr = UnmarshalDynamically<BodyTuple>::execute(rb, Ellipsoid::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   EllipsoidID e2 = static_cast<EllipsoidID>(ePtr.get());

   WALBERLA_CHECK_FLOAT_EQUAL(e1.getPosition(), e2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(e1.getLinearVel(), e2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(e1.getAngularVel(), e2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(e1.getSemiAxes(), e2->getSemiAxes());
   WALBERLA_CHECK_EQUAL(e1.getID(), e2->getID());
   WALBERLA_CHECK_EQUAL(e1.getSystemID(), e2->getSystemID());
}

void testUnion()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testUnion ***");
   UnionT u1(159, 423, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), false, false, false);
   SphereID s11 = createSphere(&u1, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), 2);
   SphereID s21 = createSphere(&u1, 4567789, Vec3(real_c(3), real_c(2), real_c(3)), real_c(1.5));
   WALBERLA_CHECK_NOT_NULLPTR( s11 );
   WALBERLA_CHECK_NOT_NULLPTR( s21 );

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, u1);
   mpi::RecvBuffer rb(sb);

   auto uPtr = UnmarshalDynamically<BodyTuple>::execute(rb, UnionT::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   UnionID u2 = static_cast<UnionID>(uPtr.get());
   WALBERLA_CHECK_NOT_NULLPTR( u2 );

   WALBERLA_CHECK_EQUAL(u1.size(), 2);
   WALBERLA_CHECK_EQUAL(u1.size(), u2->size());
   // More exhaustive tests (with inertia, rotation, etc.)
   // can be found in tests/mesh/MeshMarshalling.cpp

   //getting spheres of second union
   SphereID s12 = static_cast<SphereID> (u2->begin().getBodyID());
   SphereID s22 = static_cast<SphereID> ((++(u2->begin())).getBodyID());
   WALBERLA_CHECK_UNEQUAL( s12, s22 );

   WALBERLA_CHECK_FLOAT_EQUAL( s11->getPosition(),    s12->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL( s11->getLinearVel(),   s12->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL( s11->getAngularVel(),  s12->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL( s11->getRadius(),      s12->getRadius());
   WALBERLA_CHECK_EQUAL(       s11->getID(),          s12->getID());
   WALBERLA_CHECK_EQUAL(       s11->getSystemID(),    s12->getSystemID());

   WALBERLA_CHECK_FLOAT_EQUAL( s21->getPosition(),    s22->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL( s21->getLinearVel(),   s22->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL( s21->getAngularVel(),  s22->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL( s21->getRadius(),      s22->getRadius());
   WALBERLA_CHECK_EQUAL(       s21->getID(),          s22->getID());
   WALBERLA_CHECK_EQUAL(       s21->getSystemID(),    s22->getSystemID());
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   SetBodyTypeIDs<BodyTuple>::execute();

   testSphere();
   testBox();
   testCapsule();
   testUnion();
   testSquirmer();
   testEllipsoid();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}