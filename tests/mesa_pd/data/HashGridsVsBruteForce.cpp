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
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/HashGrids.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>

namespace walberla {
namespace mesa_pd {

using IdxPair_T = std::pair<size_t, size_t>;

// unique sorting of index pairs
bool compPair(const IdxPair_T& a, const IdxPair_T& b)
{
   if (a.first == b.first) return a.second < b.second;
   return a.first < b.first;
}

// check contact between two particles with position and radius
// mimics sphere-sphere contact detection
bool areOverlapping(Vec3 pos1, real_t radius1, Vec3 pos2, real_t radius2)
{
   return (pos2-pos1).length() < radius1 + radius2;
}

void checkTestScenario( real_t radiusRatio )
{
   const real_t radius  = real_t(0.5);
   const real_t radiusMax = std::sqrt(radiusRatio) * radius;
   const real_t radiusMin = radius / std::sqrt(radiusRatio);

   math::seedRandomGenerator( numeric_cast<std::mt19937::result_type>( 42 * walberla::mpi::MPIManager::instance()->rank() ) );

   //logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
   //logging::Logging::instance()->includeLoggingToFile("MESA_PD_Data_HashGridsVsBruteForce");
   //logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,30,30,30), // simulation domain
                                                 Vector3<uint_t>(3,3,3), // blocks in each direction
                                                 Vector3<bool>(true, true, true) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   WALBERLA_CHECK_EQUAL(forest->size(), 1);
   auto& blk = *forest->begin();

   WALBERLA_CHECK(blk.getAABB().xSize() > radiusMax && blk.getAABB().ySize() > radiusMax &&  blk.getAABB().zSize() > radiusMax,
                  "Condition for next neighbor sync violated!" );

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   data::HashGrids hg;
   std::vector<IdxPair_T> csBF;
   std::vector<IdxPair_T> csHG1;
   std::vector<IdxPair_T> csHG2;
   std::vector<IdxPair_T> csHG3;

   data::ParticleAccessor accessor(ps);

   //initialize particles
   for (int i = 0; i < 1000; ++i)
   {
      data::Particle&& p = *ps->create();
      p.getPositionRef() = Vec3( math::realRandom(blk.getAABB().xMin(), blk.getAABB().xMax()),
                                 math::realRandom(blk.getAABB().yMin(), blk.getAABB().yMax()),
                                 math::realRandom(blk.getAABB().zMin(), blk.getAABB().zMax()) );

      p.getOwnerRef() = walberla::mpi::MPIManager::instance()->rank();
      p.getInteractionRadiusRef() = math::realRandom(radiusMin, radiusMax);
   }

   //init kernels
   mpi::SyncNextNeighbors SNN;

   SNN(*ps, domain);

   ps->forEachParticlePairHalf(false,
                               kernel::SelectAll(),
                               accessor,
                               [&csBF](const size_t idx1, const size_t idx2, auto& ac)
   {
      if (areOverlapping(ac.getPosition(idx1), ac.getInteractionRadius(idx1),
                         ac.getPosition(idx2), ac.getInteractionRadius(idx2) ))
      {
         if(ac.getUid(idx2) < ac.getUid(idx1))
            csBF.push_back(IdxPair_T(idx2, idx1));
         else
            csBF.push_back(IdxPair_T(idx1, idx2));
      }
   },
   accessor );

   // insert into hash grids initially

   ps->forEachParticle(true, kernel::SelectAll(), accessor, hg, accessor);
   hg.forEachParticlePairHalf(false,
                              kernel::SelectAll(),
                              accessor,
                              [&csHG1](const size_t idx1, const size_t idx2, auto& ac)
   {
      if (areOverlapping(ac.getPosition(idx1), ac.getInteractionRadius(idx1),
                         ac.getPosition(idx2), ac.getInteractionRadius(idx2) ))
      {
         if(ac.getUid(idx2) < ac.getUid(idx1))
            csHG1.push_back(IdxPair_T(idx2, idx1));
         else
            csHG1.push_back(IdxPair_T(idx1, idx2));
      }
   },
   accessor );

   WALBERLA_CHECK_EQUAL(csBF.size(), csHG1.size());
   WALBERLA_LOG_DEVEL(csBF.size() << " contacts detected");

   std::sort(csBF.begin(), csBF.end(), compPair);
   std::sort(csHG1.begin(), csHG1.end(), compPair);

   for (size_t i = 0; i < csBF.size(); ++i)
   {
      WALBERLA_CHECK_EQUAL(csBF[i].first, csHG1[i].first);
      WALBERLA_CHECK_EQUAL(csBF[i].second, csHG1[i].second);
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("Initial insertion checked");
   WALBERLA_MPI_BARRIER();

   // redo to check clear
   hg.clear();
   ps->forEachParticle(true, kernel::SelectAll(), accessor, hg, accessor);
   hg.forEachParticlePairHalf(false,
                              kernel::SelectAll(),
                              accessor,
                              [&csHG2](const size_t idx1, const size_t idx2, auto& ac)
                              {
                                 if (areOverlapping(ac.getPosition(idx1), ac.getInteractionRadius(idx1),
                                                    ac.getPosition(idx2), ac.getInteractionRadius(idx2) ))
                                 {
                                    if(ac.getUid(idx2) < ac.getUid(idx1))
                                       csHG2.push_back(IdxPair_T(idx2, idx1));
                                    else
                                       csHG2.push_back(IdxPair_T(idx1, idx2));
                                 }
                              },
                              accessor );

   WALBERLA_CHECK_EQUAL(csBF.size(), csHG2.size());

   std::sort(csHG2.begin(), csHG2.end(), compPair);

   for (size_t i = 0; i < csBF.size(); ++i)
   {
      WALBERLA_CHECK_EQUAL(csBF[i].first, csHG2[i].first);
      WALBERLA_CHECK_EQUAL(csBF[i].second, csHG2[i].second);
   }
   WALBERLA_LOG_DEVEL_ON_ROOT("Insertion after clear checked");
   WALBERLA_MPI_BARRIER();

   // redo to check clearAll
   hg.clearAll();
   ps->forEachParticle(true, kernel::SelectAll(), accessor, hg, accessor);
   hg.forEachParticlePairHalf(false,
                              kernel::SelectAll(),
                              accessor,
                              [&csHG3](const size_t idx1, const size_t idx2, auto& ac)
                              {
                                 if (areOverlapping(ac.getPosition(idx1), ac.getInteractionRadius(idx1),
                                                    ac.getPosition(idx2), ac.getInteractionRadius(idx2) ))
                                 {
                                    if(ac.getUid(idx2) < ac.getUid(idx1))
                                       csHG3.push_back(IdxPair_T(idx2, idx1));
                                    else
                                       csHG3.push_back(IdxPair_T(idx1, idx2));
                                 }
                              },
                              accessor );

   WALBERLA_CHECK_EQUAL(csBF.size(), csHG3.size());

   std::sort(csHG3.begin(), csHG3.end(), compPair);

   for (size_t i = 0; i < csBF.size(); ++i)
   {
      WALBERLA_CHECK_EQUAL(csBF[i].first, csHG3[i].first);
      WALBERLA_CHECK_EQUAL(csBF[i].second, csHG3[i].second);
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("Insertion after clear all checked");
   WALBERLA_MPI_BARRIER();


}

/*
 * Generates particles randomly inside the domain.
 * Then checks if the Hash Grids find the same interaction pairs as the naive all-against-all check.
 * Similar to test in "LinkedCellsVsBruteForce.cpp"
 */
int main( int argc, char ** argv ) {
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_DEVEL_ON_ROOT("Checking monodisperse case");
   checkTestScenario(walberla::real_t(1.01)); // monodisperse

   WALBERLA_LOG_DEVEL_ON_ROOT("Checking polydisperse case");
   checkTestScenario(walberla::real_t(10)); // polydisperse

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
