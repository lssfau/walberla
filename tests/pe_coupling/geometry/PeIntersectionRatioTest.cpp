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
//! \file PeIntersectionRatioTest.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "pe/rigidbody/Ellipsoid.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Materials.h"

#include "pe_coupling/geometry/PeIntersectionRatio.h"

namespace pe_intersection_ratio_test
{

///////////
// USING //
///////////

using namespace walberla;

typedef std::tuple<pe::Sphere, pe::Plane, pe::Ellipsoid> BodyTypeTuple;

/*!\brief TODO
 */
//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   pe::SetBodyTypeIDs<BodyTypeTuple>::execute(); //important to be able to compare static body types in intersection function!

   const real_t epsilon( real_t(1e-5) );

   walberla::id_t sid = 0;
   walberla::id_t uid = 0;
   
   Vector3<real_t> rotationAngles( real_t(0));
   Quaternion<real_t> quat( rotationAngles );
   pe::MaterialID material = pe::Material::find("iron");


   ////////////
   // SPHERE //
   ////////////
   {
      Vector3<real_t> bodyPos(real_t(1), real_t(0), real_t(0));
      real_t radius = real_t(1);

      pe::Sphere sphere(++sid, ++uid, bodyPos, quat, radius, material, false, false, false);

      pe::RigidBody & rb = sphere; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(real_t(-0.5), real_t(0), real_t(0));
      Vector3<real_t> dir1(real_t(1), real_t(0), real_t(0));
      real_t delta1 = walberla::lbm::intersectionRatio(rb, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio with sphere wrong!");

      Vector3<real_t> pos2(real_t(1), real_t(1), real_t(1));
      Vector3<real_t> dir2(real_t(0), -real_t(1), -real_t(1));
      real_t delta2 = walberla::lbm::intersectionRatio(rb, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2) - real_t(1)) / std::sqrt(2), "Intersection ratio with sphere wrong!");
   }

   ///////////
   // PLANE //
   ///////////
   {
      Vector3<real_t> bodyPos(real_t(1), real_t(0), real_t(0));
      Vector3<real_t> bodyNormal(real_t(0), real_t(1), real_t(1));

      bodyNormal = bodyNormal.getNormalized();

      pe::Plane plane(++sid, ++uid, bodyPos, bodyNormal, bodyPos * bodyNormal, material);

      pe::RigidBody & rb = plane; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(real_t(1), real_t(0.5), real_t(0.5));
      Vector3<real_t> dir1(real_t(0), -real_t(1), -real_t(1));
      real_t delta1 = walberla::lbm::intersectionRatio(rb, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio with plane wrong!");

      Vector3<real_t> dir2(real_t(0), real_t(0), -real_t(2));
      real_t delta2 = walberla::lbm::intersectionRatio(rb, pos1, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, real_t(0.5), "Intersection ratio with plane wrong!");

      Vector3<real_t> dir3(real_t(0), -real_t(3), real_t(0));
      real_t delta3 = walberla::lbm::intersectionRatio(rb, pos1, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, real_t(1)/real_t(3), "Intersection ratio with plane wrong!");
   }

   ///////////////
   // ELLIPSOID //
   ///////////////
   {
      Vector3<real_t> bodyPos(real_t(1), real_t(0), real_t(0));
      Vector3<real_t> semiAxes1(real_t(1), real_t(1), real_t(1));

      pe::Ellipsoid ellip1(++sid, ++uid, bodyPos, quat, semiAxes1, material, false, false, false);

      pe::RigidBody & rb1 = ellip1; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(real_t(-0.5), real_t(0), real_t(0));
      Vector3<real_t> dir1(real_t(1), real_t(0), real_t(0));
      real_t delta1 = walberla::lbm::intersectionRatio(rb1, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio with ellipsoid wrong!");

      Vector3<real_t> pos2(real_t(1), real_t(1), real_t(1));
      Vector3<real_t> dir2(real_t(0), -real_t(1), -real_t(1));
      real_t delta2 = walberla::lbm::intersectionRatio(rb1, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2) - real_t(1)) / std::sqrt(2), "Intersection ratio with ellipsoid wrong!");

      Vector3<real_t> semiAxes2(real_t(2), real_t(0.5), real_t(2));
      pe::Ellipsoid ellip2(++sid, ++uid, bodyPos, quat, semiAxes2, material, false, false, false);

      pe::RigidBody & rb2 = ellip2; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos3(real_t(1), real_t(1), real_t(0));
      Vector3<real_t> dir3(real_t(0), real_t(-1), real_t(0));
      real_t delta3 = walberla::lbm::intersectionRatio(rb2, pos3, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, real_t(0.5), "Intersection ratio with ellipsoid wrong!");

   }

   return 0;

}

} //namespace pe_intersection_ratio_test

int main( int argc, char **argv ){
   pe_intersection_ratio_test::main(argc, argv);
}
