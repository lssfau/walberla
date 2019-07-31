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
//! \file   DropTestGeneral.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/GeneralContactDetection.h>

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEulerWithShape.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>

#include <core/Abort.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/waLBerlaBuildInfo.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

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

int main( int argc, char ** argv )
{
   using namespace walberla::timing;

   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<walberla::real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);

   auto p                       = ps->create();
   p->getPositionRef()          = Vec3(real_t(0), real_t(0), real_t(0));
   p->getInteractionRadiusRef() = std::numeric_limits<real_t>::infinity();
   p->getShapeIDRef()           = ss->create<data::HalfSpace>( Vec3(real_t(0), real_t(0), real_t(1)) );
   p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   p->getTypeRef()              = 0;
   using namespace walberla::mesa_pd::data::particle_flags;
   set(p->getFlagsRef(), INFINITE);
   set(p->getFlagsRef(), FIXED);
   set(p->getFlagsRef(), GLOBAL);

   auto sp                       = ps->create();
   sp->getPositionRef()          = Vec3(real_t(0), real_t(0), real_t(0.01));
   sp->getInteractionRadiusRef() = real_t(1);
   sp->getShapeIDRef()           = ss->create<data::Sphere>( real_t(0.004) );
   ss->shapes[sp->getShapeID()]->updateMassAndInertia(real_t(2707));
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[sp->getShapeID()]->getInvMass());
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[sp->getShapeID()]->getInvInertiaBF());
   sp->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   sp->getTypeRef()              = 0;

   auto bx                       = ps->create();
   bx->getPositionRef()          = Vec3(real_t(1), real_t(0), real_t(0.01));
   bx->getInteractionRadiusRef() = real_t(1);
   bx->getShapeIDRef()           = ss->create<data::Box>( Vec3(real_t(0.008*0.8)) );
   ss->shapes[bx->getShapeID()]->updateMassAndInertia(real_t(2707));
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[bx->getShapeID()]->getInvMass());
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[bx->getShapeID()]->getInvInertiaBF());
   bx->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   bx->getTypeRef()              = 0;

   auto el                       = ps->create();
   el->getPositionRef()          = Vec3(real_t(2), real_t(0), real_t(0.01));
   el->getInteractionRadiusRef() = real_t(1);
   el->getShapeIDRef()           = ss->create<data::Ellipsoid>( Vec3(real_t(0.004)) );
   ss->shapes[el->getShapeID()]->updateMassAndInertia(real_t(2707));
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[el->getShapeID()]->getInvMass());
   WALBERLA_LOG_DEVEL_VAR(ss->shapes[el->getShapeID()]->getInvInertiaBF());
   el->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   el->getTypeRef()              = 0;

   int64_t simulationSteps = 200000;
   real_t dt = real_t(0.00001);

   // Init kernels
   kernel::ExplicitEulerWithShape        explicitEulerWithShape( dt );
   kernel::SpringDashpot                 dem(1);
   auto meff = real_t(1.0) / ss->shapes[sp->getShapeID()]->getInvMass();
   dem.setParametersFromCOR(0,0,real_t(0.9), dt * real_t(20), meff);
   dem.setDampingT (0, 0, real_t(6.86e1));
   dem.setFriction (0, 0, real_t(1.2));
   WALBERLA_LOG_DEVEL_VAR(dem.getStiffness(0,0));
   WALBERLA_LOG_DEVEL_VAR(dem.getDampingN(0,0));
   WALBERLA_LOG_DEVEL_VAR(dem.calcCollisionTime(0,0,meff));
   WALBERLA_LOG_DEVEL_VAR(dem.calcCoefficientOfRestitution(0,0,meff));

   collision_detection::GeneralContactDetection gcd;
   kernel::DoubleCast                 double_cast;

   for (int64_t i=0; i < simulationSteps; ++i)
   {
      ps->forEachParticle(false,
                          kernel::SelectLocal(),
                          accessor,
                          [&](const size_t idx, auto& ac)
      {
         ac.setForce(idx,
                     Vec3(real_t(0), real_t(0), real_t(-9.81)) *
                     real_t(1.0) / ss->shapes[ac.getShapeID(idx)]->getInvMass() );
      }, accessor);

      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  accessor,
                                  [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         if (double_cast(idx1, idx2, ac, gcd, ac ))
         {
            if (ss->shapes[ac.getShapeID(idx1)]->getShapeType() == data::HalfSpace::SHAPE_TYPE)
            {
               meff = real_t(1.0) / ss->shapes[ac.getShapeID(idx2)]->getInvMass();
            } else
            {
               meff = real_t(1.0) / ss->shapes[ac.getShapeID(idx1)]->getInvMass();
            }

            dem.setParametersFromCOR(0,0,real_t(0.9), dt * real_t(20), meff);
            dem(gcd.getIdx1(), gcd.getIdx2(), ac, gcd.getContactPoint(), gcd.getContactNormal(), gcd.getPenetrationDepth());
         }
      },
      accessor );

      ps->forEachParticle(false,
                          kernel::SelectLocal(),
                          accessor,
                          explicitEulerWithShape,
                          accessor);

//      if(i%1 == 0)
//         WALBERLA_LOG_DEVEL_VAR(sp->getPosition()[2]);
   }

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sp->getPosition()[2], real_t(0.004),     real_t(1e-7));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(bx->getPosition()[2], real_t(0.004*0.8), real_t(1e-7));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(el->getPosition()[2], real_t(0.004),     real_t(1e-7));

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
