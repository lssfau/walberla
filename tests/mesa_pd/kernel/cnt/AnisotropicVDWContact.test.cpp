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

#include "SphericalSegmentAccessor.h"

#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/cnt/AnisotropicVDWContact.h"
#include "mesa_pd/kernel/cnt/ViscousDamping.h"
#include "mesa_pd/kernel/cnt/Parameters.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "vtk/VTKOutput.h"

namespace walberla {
using namespace walberla::mesa_pd;

int main(int argc, char **argv)
{
   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<walberla::real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT("loading configuration parameters");
   constexpr auto numSimulationSteps = 20000ll;
   constexpr auto outputInterval = 100ll;

   WALBERLA_LOG_INFO_ON_ROOT("creating initial particle setup");
   auto ps = std::make_shared<data::ParticleStorage>(10);
   auto ac = SphericalSegmentAccessor(ps);

   using namespace kernel::cnt;
   data::Particle &&sp1 = *ps->create();
   sp1.setPosition(Vec3(0_r, 0_r, 0_r));
   sp1.setSegmentID(1);
   sp1.setClusterID(1);

   data::Particle &&sp2 = *ps->create();
   sp2.setPosition(Vec3(20_r, 20_r, 20_r));
   sp2.setSegmentID(2);
   sp2.setClusterID(2);

   WALBERLA_LOG_INFO_ON_ROOT("setting up VTK output");
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   vtkOutput->addOutput<data::SelectParticlePosition>("position");
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput,
                                                                   "cnt",
                                                                   1,
                                                                   "vtk_integrated",
                                                                   "particles",
                                                                   false,
                                                                   false);

   WALBERLA_LOG_INFO_ON_ROOT("setting up interaction models");
   kernel::cnt::AnisotropicVDWContact vdW_anisotropic;
   kernel::cnt::ViscousDamping viscous_damping(0.1_r * 1052.0_r, 0.1_r * 1052.0_r);
   kernel::VelocityVerletPreForceUpdate vv_pre(kernel::cnt::dT);
   kernel::VelocityVerletPostForceUpdate vv_post(kernel::cnt::dT);

   WALBERLA_LOG_INFO_ON_ROOT("running simulation");

   real_t U = 0_r;
   for (auto i = 0; i < numSimulationSteps; ++i)
   {
      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_pre,
                          ac);

      U = 0_r;
      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  ac,
                                  [&](size_t p_idx1, size_t p_idx2)
                                  {
                                     vdW_anisotropic(p_idx1, p_idx2, ac);
                                     U += vdW_anisotropic.getLastEnergy();
                                     viscous_damping(p_idx1, p_idx2, ac);
                                  });

      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_post,
                          ac);

      if( i % outputInterval == 0 )
      {
//         vtkWriter->write();
//         WALBERLA_LOG_DEVEL(i << " : " << U);
      }
   }

   WALBERLA_CHECK_LESS(U, -1.16_r)

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char *argv[])
{
   return walberla::main(argc, argv);
}
