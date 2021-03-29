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

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>

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

class comp
{
public:
   comp(std::vector<collision_detection::AnalyticContactDetection>& cs) : cs_(cs) {}
   bool operator()(const size_t& c1, const size_t& c2)
   {
      if (cs_[c1].getIdx1() == cs_[c2].getIdx1()) return cs_[c1].getIdx2() < cs_[c2].getIdx2();
      return cs_[c1].getIdx1() < cs_[c2].getIdx1();
   }
   std::vector<collision_detection::AnalyticContactDetection>& cs_;
};


int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   math::seedRandomGenerator( numeric_cast<std::mt19937::result_type>( 42 * walberla::mpi::MPIManager::instance()->rank() ) );

   //logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
   //logging::Logging::instance()->includeLoggingToFile("MESA_PD_Kernel_SyncNextNeighbor");
   //logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,30,30,30), // simulation domain
                                                 Vector3<uint_t>(3,3,3), // blocks in each direction
                                                 Vector3<bool>(true, true, true) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   WALBERLA_CHECK_EQUAL(forest->size(), 1);
   const Block& blk = *static_cast<blockforest::Block*>(&*forest->begin());

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   data::LinkedCells     lc(blk.getAABB(), real_t(1.1));
   std::vector<collision_detection::AnalyticContactDetection> cs1(100);
   std::vector<collision_detection::AnalyticContactDetection> cs2(100);

   ParticleAccessorWithShape accessor(ps, ss);

   //initialize particles
   const real_t radius  = real_t(0.5);
   auto smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   for (int i = 0; i < 1000; ++i)
   {
      data::Particle&& p          = *ps->create();
      p.getPositionRef()          = Vec3( math::realRandom(blk.getAABB().xMin(), blk.getAABB().xMax()),
                                       math::realRandom(blk.getAABB().yMin(), blk.getAABB().yMax()),
                                       math::realRandom(blk.getAABB().zMin(), blk.getAABB().zMax()) );
      p.getInteractionRadiusRef() = radius;
      p.getShapeIDRef()           = smallSphere;
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   }

   //init kernels
   kernel::InsertParticleIntoLinkedCells  ipilc;
   mpi::SyncNextNeighbors                 SNN;

   SNN(*ps, domain);

   lc.clear();
   ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);

   ps->forEachParticlePairHalf(false,
                               kernel::SelectAll(),
                               accessor,
                               [&cs1](const size_t idx1, const size_t idx2, auto& ac)
   {
      collision_detection::AnalyticContactDetection         acd;
      kernel::DoubleCast               double_cast;
      if (double_cast(idx1, idx2, ac, acd, ac ))
      {
         cs1.push_back(acd);
      }
   },
   accessor );

   lc.forEachParticlePairHalf(false,
                              kernel::SelectAll(),
                              accessor,
                              [&cs2](const size_t idx1, const size_t idx2, auto& ac)
   {
      collision_detection::AnalyticContactDetection         acd;
      kernel::DoubleCast               double_cast;
      if (double_cast(idx1, idx2, ac, acd, ac ))
      {
         cs2.push_back(acd);
      }
   },
   accessor );

   WALBERLA_CHECK_EQUAL(cs1.size(), cs2.size());
   WALBERLA_LOG_DEVEL(cs1.size() << " contacts detected");

   std::vector<size_t> cs1_idx(cs1.size());
   std::vector<size_t> cs2_idx(cs2.size());
   std::iota(cs1_idx.begin(), cs1_idx.end(), 0);
   std::iota(cs2_idx.begin(), cs2_idx.end(), 0);
   std::sort(cs1_idx.begin(), cs1_idx.end(), comp(cs1));
   std::sort(cs2_idx.begin(), cs2_idx.end(), comp(cs2));


   for (size_t i = 0; i < cs1.size(); ++i)
   {
      WALBERLA_CHECK_EQUAL(cs1[cs1_idx[i]].getIdx1(), cs2[cs2_idx[i]].getIdx1());
      WALBERLA_CHECK_EQUAL(cs1[cs1_idx[i]].getIdx2(), cs2[cs2_idx[i]].getIdx2());
   }

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
