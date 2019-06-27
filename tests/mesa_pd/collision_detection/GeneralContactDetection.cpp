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
//! \file   GeneralContactDetection.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/GeneralContactDetection.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/kernel/DoubleCast.h>

#include <core/Abort.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/waLBerlaBuildInfo.h>

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

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};

void generalContactDetection()
{
   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);

   auto e0 = ps->create();
   e0->setPosition(Vec3(real_t(0),real_t(0),real_t(0)));
   e0->setShapeID(ss->create<data::Ellipsoid>(Vec3(real_t(1),real_t(2),real_t(3))));

   auto e1 = ps->create();
   e1->setPosition(Vec3(real_t(1.9),real_t(0),real_t(0)));
   e1->setShapeID(ss->create<data::Ellipsoid>(Vec3(real_t(1),real_t(2),real_t(3))));

   auto p1 = ps->create();
   p1->setPosition(Vec3(real_t(-0.9),real_t(0),real_t(0)));
   p1->setShapeID(ss->create<data::HalfSpace>(Vec3(real_t(1),real_t(0),real_t(0))));

   auto cb1 = ps->create();
   cb1->setPosition(Vec3(real_t(0),real_t(0),real_t(0)));
   cb1->setShapeID(ss->create<data::CylindricalBoundary>(real_t(3), Vec3(real_t(0),real_t(0),real_t(1))));

   collision_detection::GeneralContactDetection gcd;
   kernel::DoubleCast double_cast;

   WALBERLA_CHECK(double_cast(0, 1, accessor, gcd, accessor));
   WALBERLA_CHECK_FLOAT_EQUAL( gcd.getContactPoint(), Vec3(real_t(0.95),real_t(0),real_t(0)) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getContactNormal(), Vec3(real_t(-1),real_t(0),real_t(0)), real_t(1e-3) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getPenetrationDepth(), real_t(-0.1), real_t(1e-3)  );
   e1->setPosition(Vec3(real_t(2.1),real_t(0),real_t(0)));
   WALBERLA_CHECK(!double_cast(0, 1, accessor, gcd, accessor));

   WALBERLA_CHECK(double_cast(0, 2, accessor, gcd, accessor));
   WALBERLA_CHECK_FLOAT_EQUAL( gcd.getContactPoint(), Vec3(real_t(-0.95),real_t(0),real_t(0)) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getContactNormal(), Vec3(real_t(-1),real_t(0),real_t(0)), real_t(1e-3) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getPenetrationDepth(), real_t(-0.1), real_t(1e-3) );
   WALBERLA_CHECK(!double_cast(1, 2, accessor, gcd, accessor));

   WALBERLA_CHECK(double_cast(1, 3, accessor, gcd, accessor));
   WALBERLA_CHECK_FLOAT_EQUAL( gcd.getContactPoint(), Vec3(real_t(3.05),real_t(0),real_t(0)) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getContactNormal(), Vec3(real_t(+1),real_t(0),real_t(0)), real_t(1e-3) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( gcd.getPenetrationDepth(), real_t(-0.1), real_t(1e-3) );
   WALBERLA_CHECK(!double_cast(0, 3, accessor, gcd, accessor));
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mesa_pd::generalContactDetection();
   return EXIT_SUCCESS;
}
