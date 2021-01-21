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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Accessor.h"
#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/cnt/Parameters.h"
#include "mesa_pd/kernel/cnt/VBondContact.h"
#include "mesa_pd/kernel/VelocityVerlet.h"

#include "core/Environment.h"
#include "core/math/Constants.h"

namespace walberla {
using namespace walberla::mesa_pd;

int main(int argc, char **argv)
{
   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT("loading configuration parameters");
   constexpr auto numSimulationSteps = 500ll;

   WALBERLA_LOG_INFO_ON_ROOT("creating initial particle setup");
   auto ps = std::make_shared<data::ParticleStorage>(10);
   auto ac = Accessor(ps);

   for (auto i = 0; i < 10; ++i)
   {
      data::Particle &&p = *ps->create();
      p.setPosition(Vec3(0_r, 0_r, real_c(i) * 13.56_r));
      p.setSegmentID(i);
      p.setClusterID(1);
      if (i == 0)
         data::particle_flags::set(p.getFlagsRef(), data::particle_flags::FIXED);
      p.getRotationRef().rotate(Vec3(0_r, 1_r, 0_r), -0.5_r * math::pi);
      p.getRotationRef().rotate(Vec3(0_r, 0_r, 1_r), 0_r);
   }
   data::Particle &&last_segment = *(ps->end() - 1);

   WALBERLA_LOG_INFO_ON_ROOT("setting up interaction models");
   kernel::cnt::VBondContact vbond;
   kernel::VelocityVerletPreForceUpdate vv_pre(kernel::cnt::dT);
   kernel::VelocityVerletPostForceUpdate vv_post(kernel::cnt::dT);

   WALBERLA_LOG_INFO_ON_ROOT("running simulation");
   auto appliedForce = Vec3(1_r, 0_r, 0_r);
   auto appliedTorque = Vec3(0_r);
   std::ofstream fout("output.txt");
   for (auto i = 0; i < numSimulationSteps; ++i)
   {
      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_pre,
                          ac);

      last_segment.setForce(appliedForce);
      last_segment.setTorque(appliedTorque);
      constexpr auto cutoff2 = kernel::cnt::outer_radius * kernel::cnt::outer_radius;

      real_t tensileEnergy = 0_r;
      real_t shearEnergy = 0_r;
      real_t bendingEnergy = 0_r;
      real_t twistingEnergy = 0_r;

      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  ac,
                                  [&](size_t p_idx1, size_t p_idx2, Accessor &ac)
                                  {
                                     if ((ac.getPosition(p_idx1) - ac.getPosition(p_idx2)).sqrLength() < cutoff2)
                                     {
                                        vbond(p_idx1, p_idx2, ac);
                                        tensileEnergy += vbond.tensileEnergy;
                                        shearEnergy += vbond.shearEnergy;
                                        bendingEnergy += vbond.bendingEnergy;
                                        twistingEnergy += vbond.twistingEnergy;
                                     }
                                  },
                                  ac);

      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_post,
                          ac);

      fout << last_segment.getPosition()[0] << " "
           << tensileEnergy << " "
           << shearEnergy << " "
           << bendingEnergy << " "
           << twistingEnergy
           << std::endl;
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char *argv[])
{
   return walberla::main(argc, argv);
}
