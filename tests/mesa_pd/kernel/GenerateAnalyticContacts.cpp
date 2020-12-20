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
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/mpi/Reduce.h>

#include <iostream>
#include <memory>

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
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,3,3,3), // simulation domain
                                                 Vector3<uint_t>(3,3,3), // blocks in each direction
                                                 Vector3<bool>(true, true, true) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);
   data::LinkedCells      lc(math::AABB(-1,-1,-1,4,4,4), real_t(1));

   //initialize particles
   const real_t radius  = real_t(0.6);
   const real_t spacing = real_t(1.0);
   auto smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   WALBERLA_CHECK_EQUAL(forest->size(), 1);
   const Block& blk = *static_cast<blockforest::Block*>(&*forest->begin());

   for (auto pt : grid_generator::SCGrid(blk.getAABB(), Vec3(spacing, spacing, spacing) * real_c(0.5), spacing))
   {
      data::Particle&& p          = *ps->create();
      p.getPositionRef()          = pt;
      p.getInteractionRadiusRef() = radius;
      p.getShapeIDRef()           = smallSphere;
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   }

   //init kernels
   kernel::InsertParticleIntoLinkedCells  insert_particle_into_linked_cells;
   mpi::ReduceProperty                    RP;
   mpi::SyncNextNeighbors                 SNN;


   SNN(*ps, domain);
   ps->forEachParticlePairHalf(false,
                               kernel::ExcludeInfiniteInfinite(),
                               accessor,
                               [&](const size_t idx1, const size_t idx2, auto& ac)
   {
      collision_detection::AnalyticContactDetection         acd;
      kernel::DoubleCast               double_cast;
      mpi::ContactFilter               contact_filter;
      if (double_cast(idx1, idx2, ac, acd, ac ))
      {
         if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
         {
            ac.getForceRef(acd.getIdx1()) += acd.getContactNormal() + Vec3(1,0,0);
            ac.getForceRef(acd.getIdx2()) -= acd.getContactNormal() - Vec3(1,0,0);
         }
      }
   },
   accessor );

   RP.operator()<ForceTorqueNotification>(*ps);

   WALBERLA_CHECK_FLOAT_EQUAL(ps->getForce(0), Vec3(6,0,0));

   //TEST WITH LINKED CELLS
   ps->forEachParticle(false, kernel::SelectAll(), accessor, [](size_t idx, auto& ac) {ac.setForce(idx, Vec3(0,0,0));}, accessor);
   lc.clear();
   ps->forEachParticle(false, kernel::SelectAll(), accessor, insert_particle_into_linked_cells, accessor, lc);

   lc.forEachParticlePairHalf(false,
                              kernel::ExcludeInfiniteInfinite(),
                              accessor,
                              [&domain](const size_t idx1, const size_t idx2, auto& ac)
   {
      collision_detection::AnalyticContactDetection         acd;
      kernel::DoubleCast               double_cast;
      mpi::ContactFilter               contact_filter;
      if (double_cast(idx1, idx2, ac, acd, ac ))
      {
         if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
         {
            ac.getForceRef(acd.getIdx1()) += acd.getContactNormal() + Vec3(1,0,0);
            ac.getForceRef(acd.getIdx2()) -= acd.getContactNormal() - Vec3(1,0,0);
         }
      }
   },
   accessor );

   RP.operator()<ForceTorqueNotification>(*ps);

   WALBERLA_CHECK_FLOAT_EQUAL(ps->getForce(0), Vec3(6,0,0));

   //ps->forEachParticle(false, [](data::ParticleStorage& ps, size_t idx) {WALBERLA_LOG_DEVEL_ON_ROOT(*ps[idx]);});
   //cs.forEachContact(false, [](data::ContactStorage& cs, size_t idx) {WALBERLA_LOG_DEVEL_ON_ROOT(*cs[idx]);});

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
