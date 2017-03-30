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
//! \ingroup pe_tests
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include "pe/rigidbody/Sphere.h"
#include "pe/communication/DynamicMarshalling.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Materials.h"

#include <boost/tuple/tuple.hpp>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::pe::communication;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   typedef boost::tuple<Sphere> BodyTuple ;
   SetBodyTypeIDs<BodyTuple>::execute();

   MaterialID iron = Material::find("iron");

   SphereID s1 = new Sphere(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), Vec3(0,0,0), Quat(), 5, iron, false, false, false);
   s1->setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   s1->setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, *s1);
   mpi::RecvBuffer rb(sb);

   SphereID s2 = static_cast<SphereID> (UnmarshalDynamically<BodyTuple>::execute(rb, Sphere::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100))));

   WALBERLA_CHECK_FLOAT_EQUAL(s1->getPosition(), s2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(s1->getLinearVel(), s2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(s1->getAngularVel(), s2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(s1->getRadius(), s2->getRadius());
   WALBERLA_CHECK_EQUAL(s1->getID(), s2->getID());
   WALBERLA_CHECK_EQUAL(s1->getSystemID(), s2->getSystemID());

   return EXIT_SUCCESS;
}
