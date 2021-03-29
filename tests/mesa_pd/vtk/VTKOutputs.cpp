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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/data/ParticleStorage.h>

#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <vtk/VTKOutput.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();
   WALBERLA_CHECK_EQUAL(mpi::MPIManager::instance()->numProcesses(), 8, "please run with 8 processes");

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage> (100);

   const math::AABB domain(real_t(0), real_t(0), real_t(0), real_t(4), real_t(4), real_t(4));
   auto forest = blockforest::createBlockForest( domain,
                                                 Vector3<uint_t>(2,2,2),
                                                 Vector3<bool>(false, false, false) );

   auto rank = mpi::MPIManager::instance()->rank();

   //initialize particles
   const Vec3 shift(real_t(0.5), real_t(0.5), real_t(0.5));
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), shift, real_t(1)))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p = ps->create();
         p->setPosition( pt );
         p->setOwner( mpi::MPIManager::instance()->rank() );
         p->getRotationRef().rotate(Vec3(0_r,1_r,0_r), real_c(rank) * math::pi * 0.5_r);
      }
   }

   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");
   vtkOutput->addOutput<data::SelectParticleLinearVelocity>("velocity");
   vtkOutput->addOutput<data::SelectParticleRotation>("rotation");
   vtkOutput->setParticleSelector( [rank](const data::ParticleStorage::iterator& pIt) {return pIt->getIdx() < uint_c(rank);} );
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "particles", 1, "vtk_outputs", "simulation_step", false, false);

   vtkWriter->write();

   WALBERLA_CHECK_EQUAL(vtkOutput->getParticlesWritten(), rank);

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
