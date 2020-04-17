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

#include "mesh/pe/rigid_body/ConvexPolyhedron.h"
#include "mesh/pe/rigid_body/ConvexPolyhedronFactory.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/QHull.h"
#include "mesh/pe/Types.h"
#include "pe/rigidbody/Squirmer.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/rigidbody/Union.h"
#include "pe/communication/rigidbody/Squirmer.h"
#include "pe/communication/DynamicMarshalling.h"
#include "mesh/pe/communication/ConvexPolyhedron.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Materials.h"

#include <memory>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::pe::communication;

using UnionT = Union<mesh::pe::ConvexPolyhedron>;
using UnionID = UnionT *;
using UnionPtr = std::unique_ptr<UnionT>;

typedef std::tuple<mesh::pe::ConvexPolyhedron, UnionT> BodyTuple ;

std::vector<Vector3<real_t>> generateOctahedron( const real_t radius)
{

   std::vector<Vector3<real_t>> okta( 6 );
   for(size_t i = 0; i < 6; i++){
      auto &p = okta[i];
      p[i%3]=(i<3) ? radius: -radius;
   }
   return okta;
}

// Checks two mesh::TriangleMesh for pseudo-equality
void checkMeshEquals(const mesh::TriangleMesh &m1, const mesh::TriangleMesh &m2){
	// Very basic checks
	WALBERLA_CHECK_FLOAT_EQUAL(mesh::computeVolume(m1), mesh::computeVolume(m2));
	WALBERLA_CHECK_FLOAT_EQUAL(mesh::toWalberla(mesh::computeCentroid(m1)), mesh::toWalberla(mesh::computeCentroid(m2)));
	WALBERLA_CHECK_FLOAT_EQUAL(mesh::computeInertiaTensor(m1), mesh::computeInertiaTensor(m2));
}

// Checks two convexPolyhedrons for equality
void checkConvexPolyhedronEquals(const mesh::pe::ConvexPolyhedron &b1, const mesh::pe::ConvexPolyhedron &b2){
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getPosition(), b2.getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getLinearVel(), b2.getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getAngularVel(), b2.getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getInertia(), b2.getInertia());
   WALBERLA_CHECK_EQUAL(b1.getMaterial(), b2.getMaterial());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getQuaternion(), b2.getQuaternion());
   WALBERLA_CHECK_FLOAT_EQUAL(b1.getRotation(), b2.getRotation());
   // Check equality of the meshes
   checkMeshEquals(b1.getMesh(), b2.getMesh());
   WALBERLA_CHECK_EQUAL(b1.getID(), b2.getID());
   WALBERLA_CHECK_EQUAL(b1.getSystemID(), b2.getSystemID());
}

void testConvexPolyhedron()
{
   WALBERLA_LOG_INFO_ON_ROOT("*** testConvexPolyhedron ***");

   // Generate mesh
   shared_ptr< mesh::TriangleMesh > octamesh = make_shared<mesh::TriangleMesh>();
   mesh::QHull< mesh::TriangleMesh > qhull( generateOctahedron(real_t(1.0)), octamesh );
   qhull.run();
   
   MaterialID iron = Material::find("iron");
   
   mesh::pe::ConvexPolyhedron b1(759846, 1234794, Vec3(real_c(1), real_c(2), real_c(3)), Quat(), *octamesh, iron, false, true, false);
   b1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   b1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));

   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, b1);
   mpi::RecvBuffer rb(sb);

   auto bPtr = UnmarshalDynamically<BodyTuple>::execute(rb, mesh::pe::ConvexPolyhedron::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   auto b2 = dynamic_cast<mesh::pe::ConvexPolyhedronID>(bPtr.get());
   checkConvexPolyhedronEquals(b1, *b2);
   
}

void testUnion()
{
   
   WALBERLA_LOG_INFO_ON_ROOT("*** testUnion ***");
   // Generate mesh
   shared_ptr< mesh::TriangleMesh > octamesh = make_shared<mesh::TriangleMesh>();
   mesh::QHull< mesh::TriangleMesh > qhull( generateOctahedron(real_t(1.0)), octamesh );
   qhull.run();
   
   MaterialID iron = Material::find("iron");
   Quat q(real_t(0.3), real_t(0.1), real_t(0.9), real_t(0.3));
   Quat w(real_t(0.1), real_t(0.3), real_t(0.9), real_t(0.3));
   UnionT u1(159, 423, Vec3(real_c(1), real_c(2), real_c(3)),  q, false, false, false);
   u1.add(std::make_unique<mesh::pe::ConvexPolyhedron>(753326, 1267824, Vec3(real_c(2), real_c(2), real_c(3)), w, *octamesh, iron, false, true, false));
   u1.add(std::make_unique<mesh::pe::ConvexPolyhedron>(753246, 1233424, Vec3(real_c(-1), real_c(4), real_c(-2)), w, *octamesh, iron, false, true, false));

   u1.setLinearVel(Vec3(real_c(5.2), real_c(6.3), real_c(7.4)));
   u1.setAngularVel(Vec3(real_c(1.2), real_c(2.3), real_c(3.4)));
   
   mpi::SendBuffer sb;
   MarshalDynamically<BodyTuple>::execute(sb, u1);
   mpi::RecvBuffer rb(sb);

   auto uPtr = UnmarshalDynamically<BodyTuple>::execute(rb, UnionT::getStaticTypeID(), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)), math::AABB(Vec3(-100,-100,-100), Vec3(100,100,100)));
   auto u2 = dynamic_cast<UnionID>(uPtr.get());
   WALBERLA_CHECK_NOT_NULLPTR( u2 );
   WALBERLA_CHECK_EQUAL(u1.size(), 2);
   WALBERLA_CHECK_EQUAL(u1.size(), u2->size());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getInertia(), u2->getInertia());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getPosition(), u2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getLinearVel(), u2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getAngularVel(), u2->getAngularVel());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getQuaternion(), u2->getQuaternion());
   WALBERLA_CHECK_FLOAT_EQUAL(u1.getRotation(), u2->getRotation());

   //getting polyhedrons of first union
   mesh::pe::ConvexPolyhedronID p11 = dynamic_cast<mesh::pe::ConvexPolyhedronID > (u1.begin().getBodyID());
   mesh::pe::ConvexPolyhedronID p21 = dynamic_cast<mesh::pe::ConvexPolyhedronID > ((++(u1.begin())).getBodyID());
   
   //getting polyhedrons of second union
   mesh::pe::ConvexPolyhedronID p12 = dynamic_cast<mesh::pe::ConvexPolyhedronID > (u2->begin().getBodyID());
   mesh::pe::ConvexPolyhedronID p22 = dynamic_cast<mesh::pe::ConvexPolyhedronID > ((++(u2->begin())).getBodyID());
   
   checkConvexPolyhedronEquals(*p11, *p12);
   checkConvexPolyhedronEquals(*p21, *p22);

}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   SetBodyTypeIDs<BodyTuple>::execute();
   testConvexPolyhedron();
   testUnion();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
