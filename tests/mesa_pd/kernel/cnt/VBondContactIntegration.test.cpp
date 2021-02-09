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
#include "mesa_pd/kernel/cnt/Parameters.h"
#include "mesa_pd/kernel/cnt/VBondContact.h"
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

   if (std::is_same<real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT("loading configuration parameters");
   constexpr auto numSimulationSteps = 500ll;

   WALBERLA_LOG_INFO_ON_ROOT("creating initial particle setup");
   auto ps = std::make_shared<data::ParticleStorage>(10);
   auto ac = SphericalSegmentAccessor(ps);

   for (auto i = 0; i < 10; ++i)
   {
      data::Particle &&p = *ps->create();
      p.setPosition(Vec3(500_r, 500_r, 500_r + real_c(i) * 13.56_r));
      p.setSegmentID(i);
      p.setClusterID(1);
      if (i == 0)
         data::particle_flags::set(p.getFlagsRef(), data::particle_flags::FIXED);
      p.getRotationRef().rotate(Vec3(0_r, 1_r, 0_r), -0.5_r * math::pi);
      p.getRotationRef().rotate(Vec3(0_r, 0_r, 1_r), 0_r);
   }
   data::Particle &&last_segment = *(ps->end() - 1);

   WALBERLA_LOG_INFO_ON_ROOT("setting up VTK output");
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   vtkOutput->addOutput<data::SelectParticlePosition>("position");
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput,
                                                                   "cnt",
                                                                   1,
                                                                   "vtk",
                                                                   "particles",
                                                                   false,
                                                                   false);

   WALBERLA_LOG_INFO_ON_ROOT("setting up interaction models");
   kernel::cnt::VBondContact vbond;
   kernel::VelocityVerletPreForceUpdate vv_pre(kernel::cnt::dT);
   kernel::VelocityVerletPostForceUpdate vv_post(kernel::cnt::dT);

   WALBERLA_LOG_INFO_ON_ROOT("running simulation");
   auto appliedForce = Vec3(1_r, 0_r, 0_r);
   auto appliedTorque = Vec3(0_r, 0_r, 1_r);

   real_t tensileEnergy = 0_r;
   real_t shearEnergy = 0_r;
   real_t bendingEnergy = 0_r;
   real_t twistingEnergy = 0_r;

   std::ofstream fout("VBondContactEnergies.txt");
   for (auto i = 0; i < numSimulationSteps; ++i)
   {
      vtkWriter->write();

      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_pre,
                          ac);

      last_segment.setForce(appliedForce);
      last_segment.setTorque(appliedTorque);
      constexpr auto cutoff2 = kernel::cnt::outer_radius * kernel::cnt::outer_radius;

      tensileEnergy = 0_r;
      shearEnergy = 0_r;
      bendingEnergy = 0_r;
      twistingEnergy = 0_r;

      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  ac,
                                  [&](size_t p_idx1, size_t p_idx2)
                                  {
                                     if (ac.getClusterID(p_idx1) != ac.getClusterID(p_idx2)) return;
                                     if (std::abs(ac.getSegmentID(p_idx1) - ac.getSegmentID(p_idx2)) != 1) return;
                                     if ((ac.getPosition(p_idx1) - ac.getPosition(p_idx2)).sqrLength() < cutoff2)
                                     {
                                        vbond(p_idx1, p_idx2, ac);
                                        tensileEnergy += vbond.getLastTensileEnergy();
                                        shearEnergy += vbond.getLastShearEnergy();
                                        bendingEnergy += vbond.getLastBendingEnergy();
                                        twistingEnergy += vbond.getLastTwistingEnergy();
                                     }
                                  });

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

   WALBERLA_CHECK_FLOAT_EQUAL(tensileEnergy,  1.88111638964774328e-02_r);
   WALBERLA_CHECK_FLOAT_EQUAL(shearEnergy,    2.04795345750102054e-01_r);
   WALBERLA_CHECK_FLOAT_EQUAL(bendingEnergy,  3.28859587360327978e+01_r);
   WALBERLA_CHECK_FLOAT_EQUAL(twistingEnergy, 1.44177931971837983e-03_r);

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char *argv[])
{
   return walberla::main(argc, argv);
}
