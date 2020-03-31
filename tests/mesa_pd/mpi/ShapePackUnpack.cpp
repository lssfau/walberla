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
//! \file   PackUnpack.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include <mesa_pd/mpi/ShapePackUnpack.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <iostream>
#include <memory>

namespace walberla {

using namespace walberla::mesa_pd;

void checkIdenticalOrInf(const Mat3& mat0, const Mat3& mat1) {
   for (uint_t i = 0; i <= 8; i++) {
      WALBERLA_CHECK((math::isinf(mat0[i]) && math::isinf(mat1[i])) || realIsIdentical(mat0[i], mat1[i]), "Matrices don't match: " << mat0 << " vs. " << mat1);
   }
}

void checkIdenticalOrInf(real_t a, real_t b) {
   WALBERLA_CHECK((math::isinf(a) && math::isinf(b)) || realIsIdentical(a, b), "Values don't match: " << a << " vs. " << b);
}

void checkBox()
{
   using namespace walberla::mpi;

   std::shared_ptr<data::BaseShape> bs0 = std::make_shared<data::Box>(Vec3(real_t(2.53), real_t(4.53), real_t(3.53)));
   bs0->updateMassAndInertia(real_t(123));
   std::shared_ptr<data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("packing box");
   SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("unpacking box");
   RecvBuffer rb(sb);
   rb >> bs1;

   WALBERLA_LOG_INFO("checking box");
   WALBERLA_CHECK_EQUAL(bs0->getShapeType(), data::Box::SHAPE_TYPE);
   WALBERLA_CHECK_EQUAL(bs1->getShapeType(), data::Box::SHAPE_TYPE);
   WALBERLA_CHECK_IDENTICAL(bs0->getMass(), bs1->getMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   checkIdenticalOrInf(bs0->getInertiaBF(), bs1->getInertiaBF());
   checkIdenticalOrInf(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto bx0 = static_cast<data::Box*> (bs0.get());
   auto bx1 = static_cast<data::Box*> (bs1.get());

   WALBERLA_CHECK_IDENTICAL(bx0->getEdgeLength(), bx1->getEdgeLength());
}

void checkCylindricalBoundary()
{
   using namespace walberla::mpi;

   std::shared_ptr<data::BaseShape> bs0 = std::make_shared<data::CylindricalBoundary>(real_t(9.99),
                                                                                      Vec3(real_t(2.53), real_t(4.53), real_t(3.53)).getNormalized());
   bs0->updateMassAndInertia(real_t(123));
   std::shared_ptr<data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("packing cylindrical boundary");
   SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("unpacking cylindrical boundary");
   RecvBuffer rb(sb);
   rb >> bs1;

   WALBERLA_LOG_INFO("checking cylindrical boundary");
   WALBERLA_CHECK_EQUAL(bs0->getShapeType(), data::CylindricalBoundary::SHAPE_TYPE);
   WALBERLA_CHECK_EQUAL(bs1->getShapeType(), data::CylindricalBoundary::SHAPE_TYPE);
   WALBERLA_CHECK_IDENTICAL(bs0->getMass(), bs1->getMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   checkIdenticalOrInf(bs0->getInertiaBF(), bs1->getInertiaBF());
   checkIdenticalOrInf(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto cb0 = static_cast<data::CylindricalBoundary*> (bs0.get());
   auto cb1 = static_cast<data::CylindricalBoundary*> (bs1.get());

   WALBERLA_CHECK_IDENTICAL(cb0->getAxis(), cb1->getAxis());
   WALBERLA_CHECK_IDENTICAL(cb0->getRadius(), cb1->getRadius());
}

void checkEllipsoid()
{
   using namespace walberla::mpi;

   std::shared_ptr<data::BaseShape> bs0 = std::make_shared<data::Ellipsoid>(Vec3(real_t(2.53), real_t(4.53), real_t(3.53)));
   bs0->updateMassAndInertia(real_t(123));
   std::shared_ptr<data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("packing ellipsoid");
   SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("unpacking ellipsoid");
   RecvBuffer rb(sb);
   rb >> bs1;

   WALBERLA_LOG_INFO("checking ellipsoid");
   WALBERLA_CHECK_EQUAL(bs0->getShapeType(), data::Ellipsoid::SHAPE_TYPE);
   WALBERLA_CHECK_EQUAL(bs1->getShapeType(), data::Ellipsoid::SHAPE_TYPE);
   WALBERLA_CHECK_IDENTICAL(bs0->getMass(), bs1->getMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   checkIdenticalOrInf(bs0->getInertiaBF(), bs1->getInertiaBF());
   checkIdenticalOrInf(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto el0 = static_cast<data::Ellipsoid*> (bs0.get());
   auto el1 = static_cast<data::Ellipsoid*> (bs1.get());

   WALBERLA_CHECK_IDENTICAL(el0->getSemiAxes(), el1->getSemiAxes());
}

void checkHalfSpace()
{
   using namespace walberla::mpi;

   std::shared_ptr<data::BaseShape> bs0 = std::make_shared<data::HalfSpace>(Vec3(real_t(2.53), real_t(4.53), real_t(3.53)));
   bs0->updateMassAndInertia(real_t(123));
   std::shared_ptr<data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("packing half space");
   SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("unpacking half space");
   RecvBuffer rb(sb);
   rb >> bs1;

   WALBERLA_LOG_INFO("checking half space");
   WALBERLA_CHECK_EQUAL(bs0->getShapeType(), data::HalfSpace::SHAPE_TYPE);
   WALBERLA_CHECK_EQUAL(bs1->getShapeType(), data::HalfSpace::SHAPE_TYPE);
   checkIdenticalOrInf(bs0->getMass(), bs1->getMass()); //INF
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   checkIdenticalOrInf(bs0->getInertiaBF(), bs1->getInertiaBF()); //INF
   WALBERLA_CHECK_IDENTICAL(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto hs0 = static_cast<data::HalfSpace*> (bs0.get());
   auto hs1 = static_cast<data::HalfSpace*> (bs1.get());

   WALBERLA_CHECK_IDENTICAL(hs0->getNormal(), hs1->getNormal());
}

void checkSphere()
{
   using namespace walberla::mpi;

   std::shared_ptr<data::BaseShape> bs0 = std::make_shared<data::Sphere>(real_t(2.53));
   bs0->updateMassAndInertia(real_t(123));
   std::shared_ptr<data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("packing sphere");
   SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("unpacking sphere");
   RecvBuffer rb(sb);
   rb >> bs1;

   WALBERLA_LOG_INFO("checking sphere");
   WALBERLA_CHECK_EQUAL(bs0->getShapeType(), data::Sphere::SHAPE_TYPE);
   WALBERLA_CHECK_EQUAL(bs1->getShapeType(), data::Sphere::SHAPE_TYPE);
   WALBERLA_CHECK_IDENTICAL(bs0->getMass(), bs1->getMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInertiaBF(), bs1->getInertiaBF());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto sp0 = static_cast<data::Sphere*> (bs0.get());
   auto sp1 = static_cast<data::Sphere*> (bs1.get());

   WALBERLA_CHECK_IDENTICAL(sp0->getRadius(), sp1->getRadius());
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   checkBox();
   checkCylindricalBoundary();
   checkEllipsoid();
   checkHalfSpace();
   checkSphere();

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}