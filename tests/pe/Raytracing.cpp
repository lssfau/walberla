#include "pe/utility/BodyCast.h"

#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <pe/raytracing/Ray.h>
#include <pe/raytracing/Intersects.h>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

typedef boost::tuple<Box, Capsule, Plane, Sphere> BodyTuple;

void SphereIntersectsTest()
{
   MaterialID iron = Material::find("iron");
   Sphere sp1(123, 1, Vec3(3,3,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);

   // ray through the center
   Ray ray1(Vec3(-5,3,3), Vec3(1,0,0));
   
   real_t t;
   
   WALBERLA_LOG_INFO("RAY -> SPHERE: through center (hitting)");
   WALBERLA_CHECK(intersects(&sp1, &ray1, &t));
   WALBERLA_CHECK(realIsEqual(t, real_t(6)))

   // ray tangential
   Ray ray2(Vec3(-5,3,3), Vec3(7.5,0,sqrt(real_t(15))/real_t(2)));
   
   WALBERLA_LOG_INFO("RAY -> SPHERE: tangential (not hitting)");
   WALBERLA_CHECK(!intersects(&sp1, &ray2, &t));
}


int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   
   SetBodyTypeIDs<BodyTuple>::execute();
   
   SphereIntersectsTest();
   
   return EXIT_SUCCESS;
}
