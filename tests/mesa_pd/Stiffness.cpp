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
//! \file   Stiffness.cpp
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithBaseShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ContactFilter.h"

namespace walberla
{

using namespace mesa_pd;

int main(int argc, char** argv)
{
   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   WALBERLA_CHECK(MPIManager::instance()->numProcesses() == 1)

   // Config
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   WALBERLA_LOG_INFO_ON_ROOT(*cfg);
   const Config::BlockHandle config = cfg->getBlock("Stiffness");

   const Vec3 domainSize_SI             = config.getParameter< Vec3 >("domainSize_SI");
   const real_t diameter_SI             = config.getParameter< real_t >("diameter_SI");
   const real_t densityParticle_SI      = config.getParameter< real_t >("densityParticle_SI");
   const real_t dt_SI                   = config.getParameter< real_t >("dt_SI");
   const uint_t timeSteps               = config.getParameter< uint_t >("timeSteps");
   const real_t force_SI                = config.getParameter< real_t >("force_SI");
   const real_t normalSpringConstant_SI = config.getParameter< real_t >("normalSpringConstant_SI");

   // BlockForest
   const math::AABB simulationDomain_SI(real_t(0.0), real_t(0.0), real_t(0.0), domainSize_SI[0], domainSize_SI[1],
                                        domainSize_SI[2]);

   shared_ptr< BlockForest > forest =
      blockforest::createBlockForest(simulationDomain_SI, Vec3(uint(1)), Vector3< bool >(false));
   auto domain = std::make_shared< mesa_pd::domain::BlockForestDomain >(forest);

   // MesaPD data structures
   auto ps = std::make_shared< data::ParticleStorage >(1);
   data::ParticleAccessorWithBaseShape accessor(ps);

   // Init sphere 0
   auto p0                       = ps->create();
   p0->getPositionRef()          = simulationDomain_SI.center() - Vec3(diameter_SI / 2, real_t(0), real_t(0));
   p0->getInteractionRadiusRef() = diameter_SI * real_t(0.5);
   p0->getBaseShapeRef()         = std::make_shared< data::Sphere >(p0->getInteractionRadius());
   p0->getBaseShapeRef()->updateMassAndInertia(densityParticle_SI);
   p0->getOwnerRef() = walberla::mpi::MPIManager::instance()->rank();
   p0->getTypeRef()  = 0;
   auto idxp0        = p0->getIdx();

   // Init sphere 1
   auto p1                       = ps->create();
   p1->getPositionRef()          = simulationDomain_SI.center() + Vec3(diameter_SI / 2, real_t(0), real_t(0));
   p1->getInteractionRadiusRef() = diameter_SI * real_t(0.5);
   p1->getBaseShapeRef()         = std::make_shared< data::Sphere >(p1->getInteractionRadius());
   p1->getBaseShapeRef()->updateMassAndInertia(densityParticle_SI);
   p1->getOwnerRef() = walberla::mpi::MPIManager::instance()->rank();
   p1->getTypeRef()  = 0;
   auto idxp1        = p1->getIdx();

   auto overlap = diameter_SI - (accessor.getPosition(idxp1)[0] - accessor.getPosition(idxp0)[0]);

   // Init kernels
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(dt_SI);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(dt_SI);
   kernel::LinearSpringDashpot dem(1);
   dem.setStiffnessN(0, 0, normalSpringConstant_SI);

   for (uint_t i = 0; i < timeSteps; ++i)
   {
      ps->forEachParticle(false, kernel::SelectLocal(), accessor, vvIntegratorPreForce, accessor);

      p0->setForce(Vec3(force_SI * real_t(i) / real_t(timeSteps), real_t(0), real_t(0)));
      p1->setForce(Vec3(-force_SI * real_t(i) / real_t(timeSteps), real_t(0), real_t(0)));

      ps->forEachParticlePairHalf(
         false, kernel::ExcludeInfiniteInfinite(), accessor,
         [domain, &dem, dt_SI](const size_t idx1, const size_t idx2, auto& ac) {
            kernel::DoubleCast double_cast;
            mesa_pd::mpi::ContactFilter contact_filter;
            collision_detection::AnalyticContactDetection acd;

            if (double_cast(idx1, idx2, ac, acd, ac))
            {
               if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *domain))
               {
                  dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                      acd.getPenetrationDepth(), dt_SI);
               }
            }
         },
         accessor);

      overlap = diameter_SI - (accessor.getPosition(idxp1)[0] - accessor.getPosition(idxp0)[0]);

      ps->forEachParticle(false, kernel::SelectLocal(), accessor, vvIntegratorPostForce, accessor);
   }

   WALBERLA_LOG_DEVEL_VAR(overlap)
   const real_t expectedOverlap = force_SI / normalSpringConstant_SI;
   WALBERLA_LOG_DEVEL_VAR(expectedOverlap)
   WALBERLA_CHECK_FLOAT_EQUAL(overlap, expectedOverlap)

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
