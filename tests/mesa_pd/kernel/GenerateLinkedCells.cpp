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
//! \file   GenerateLinkedCells.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/kernel/ExplicitEuler.h>
#include <mesa_pd/kernel/ForceLJ.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>

#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

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
   domain.scale(real_t(10));

   //init data structures
   auto storage = std::make_shared<data::ParticleStorage>(100);
   data::LinkedCells     linkedCells(domain.getScaled(2), real_t(1));

   data::ParticleAccessor accessor(storage);

   //initialize particles
   for (auto it = grid_generator::SCIterator(domain, Vec3(spacing, spacing, spacing) * real_c(0.5), spacing);
        it != grid_generator::SCIterator();
        ++it)
   {
      data::Particle&& p    = *storage->create();
      p.getPositionRef()    = (*it);
   }

   //init kernels
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::ForceLJ lj(1);
   kernel::ExplicitEuler integrator( real_t(0.01) );

   //timeloop
   for (auto timestep = 0; timestep < 100; ++timestep)
   {
      linkedCells.clear();
      storage->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, linkedCells);

      int particleCounter = 0;
      for (int x = 0; x < linkedCells.numCellsPerDim_[0]; ++x)
         for (int y = 0; y < linkedCells.numCellsPerDim_[1]; ++y)
            for (int z = 0; z < linkedCells.numCellsPerDim_[2]; ++z)
            {
               const uint_t cell_idx = getCellIdx(linkedCells, x, y, z);
               auto aabb = getCellAABB(linkedCells, x, y, z);
               int p_idx = linkedCells.cells_[cell_idx];
               while (p_idx != -1)
               {
                  ++particleCounter;
                  WALBERLA_CHECK( aabb.contains( storage->getPosition(uint_c(p_idx)) ),
                                  "Particle(" << p_idx << ") with position (" <<
                                  storage->getPosition(uint_c(p_idx)) <<
                                  ") not contained in cell(" << x << ", " << y << ", " << z <<
                                  ") with aabb: " << aabb << ".");
                  p_idx = storage->getNextParticle(uint_c(p_idx));
               }
            }
      WALBERLA_CHECK_EQUAL(particleCounter, storage->size());
      linkedCells.forEachParticlePairHalf(true, kernel::SelectAll(), accessor, lj, accessor);
      std::for_each(storage->begin(), storage->end(), [](data::Particle&& p){ p.getForceRef() = -p.getPosition() + p.getForce(); });
      storage->forEachParticle(true, kernel::SelectAll(), accessor, integrator, accessor);
   }

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
