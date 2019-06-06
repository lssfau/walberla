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
//! \file   01_LennardJones.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEuler.h>
#include <mesa_pd/kernel/ForceLJ.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/VelocityVerlet.h>

#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
#include <vtk/VTKOutput.h>

#include <iostream>
#include <memory>

using namespace walberla;
using namespace walberla::mesa_pd;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //domain setup
   const real_t spacing = real_t(1.0);
   math::AABB domain( Vec3(real_t(-0.5), real_t(-0.5), real_t(-0.5)),
                      Vec3(real_t(+0.5), real_t(+0.5), real_t(+0.5)));
   domain.scale(real_t(20));

   //init data structures
   auto storage = std::make_shared<data::ParticleStorage> (100);
   data::LinkedCells     linkedCells(domain.getScaled(2), real_t(1));
   data::ParticleAccessor ac(storage);

   //initialize particles
   for (auto it = grid_generator::SCIterator(domain, Vec3(spacing, spacing, spacing) * real_c(0.5), spacing);
        it != grid_generator::SCIterator();
        ++it)
   {
      data::Particle&& p = *storage->create();
      p.getPositionRef() = (*it);

      p.getLinearVelocityRef() = Vec3(math::realRandom(real_t(-1), real_t(+1)),
                                      math::realRandom(real_t(-1), real_t(+1)),
                                      math::realRandom(real_t(-1), real_t(+1)));
   }

   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(storage) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk", "simulation_step", false, false);

   // Init kernels
   kernel::ForceLJ                        lj(1);
   lj.setSigma(0,0,real_t(1));
   lj.setEpsilon(0,0,real_t(1));
   kernel::InsertParticleIntoLinkedCells  ipilc;
   kernel::VelocityVerletPreForceUpdate   vvPreForce( real_t(0.01) );
   kernel::VelocityVerletPostForceUpdate  vvPostForce( real_t(0.01) );

   for (auto timestep = 0; timestep < 1000; ++timestep)
   {
      WALBERLA_LOG_DEVEL(timestep);
      linkedCells.clear();
      storage->forEachParticle(true, kernel::SelectAll(), ac, ipilc, ac, linkedCells);
      storage->forEachParticle(true, kernel::SelectLocal(), ac, vvPreForce, ac);
      linkedCells.forEachParticlePairHalf(true, kernel::SelectAll(), ac, lj, ac);
      const real_t coeff = real_t(0.2);
      storage->forEachParticle(true,
                               kernel::SelectLocal(),
                               ac,
                               [coeff](const size_t idx, auto& access){ access.setForce(idx, -coeff*access.getPosition(idx) + access.getForce(idx)); },
                               ac);
      storage->forEachParticle(true, kernel::SelectLocal(), ac, vvPostForce, ac);
      vtkWriter->write();
   }

   return EXIT_SUCCESS;
}
