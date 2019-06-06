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
//! \file   ReduceContactHistory.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceContactHistory.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <mesa_pd/mpi/notifications/ContactHistoryNotification.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {
namespace mesa_pd {

const real_t radius = real_t(0.7);
const real_t spacing = real_t(1.0);

class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps)
      : ParticleAccessor(ps)
   {}

   data::BaseShape* getShape(const size_t /*p_idx*/) {return &sp;}
private:
   data::Sphere sp = data::Sphere(radius);
};

class IncreaseContactHistory
{
public:
   IncreaseContactHistory() = default;

   template <typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor& ac,
                   const Vec3& /*contactPoint*/,
                   const Vec3& /*contactNormal*/,
                   const real_t& /*penetrationDepth*/) const
   {
      auto& oldCH1 = ac.getOldContactHistoryRef(p_idx1)[ac.getUid(p_idx2)];
      auto& oldCH2 = ac.getOldContactHistoryRef(p_idx2)[ac.getUid(p_idx1)];
      auto& newCH1 = ac.getNewContactHistoryRef(p_idx1)[ac.getUid(p_idx2)];
      auto& newCH2 = ac.getNewContactHistoryRef(p_idx2)[ac.getUid(p_idx1)];
      newCH1.setTangentialSpringDisplacement(oldCH1.getTangentialSpringDisplacement() + Vec3(ac.getUid(p_idx2)));
      newCH2.setTangentialSpringDisplacement(oldCH2.getTangentialSpringDisplacement() + Vec3(ac.getUid(p_idx1)));
   }

};

void createSphere(data::ParticleStorage& ps, domain::IDomain& domain, const Vec3& pos)
{
   static uint64_t uid = 0;

   auto owned = domain.isContainedInProcessSubdomain( uint_c(walberla::mpi::MPIManager::instance()->rank()), pos );
   if (owned)
   {
      data::Particle&& p          = *ps.create(uid); // DO NOT DO THIS IF YOU DONT KNOW WHAT YOU ARE DOING
      p.getPositionRef()          = pos;
      p.getInteractionRadiusRef() = radius;
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
      p.getTypeRef()              = 0;
   }
   ++uid;
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,4,4,4), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(true, true, true) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   ParticleAccessorWithShape accessor(ps);

   for (auto pt : grid_generator::SCGrid(forest->getDomain(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.4), spacing))
   {
      createSphere(*ps, domain, pt);
   }

   IncreaseContactHistory                INC;
   mpi::ReduceContactHistory             ReduceContactHistory;
   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighbors                SNN;

   SNN(*ps, domain);

   const int simulationSteps = 10;
   for (int i=0; i < simulationSteps; ++i)
   {
      //if (i % visSpacing == 0)
      //{
      //   vtkWriter->write();
      //}


      ps->forEachParticlePairHalf(true,
                                  kernel::SelectAll(),
                                  accessor,
                                  [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         collision_detection::AnalyticContactDetection acd;
         kernel::DoubleCast double_cast;
         const mpi::ContactFilter contact_filter;
         if (double_cast(idx1, idx2, ac, acd, ac ))
         {
            if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
            {
               INC(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
            }
         }
      },
      accessor );

      ReduceContactHistory(*ps);

      SNN(*ps, domain);
   }


   for (data::Particle&& p : *ps)
   {
      WALBERLA_CHECK_EQUAL(p.getOldContactHistory().size(), 6, p);
      WALBERLA_CHECK_EQUAL(p.getNewContactHistory().size(), 0, p);
      for (auto& entry : p.getOldContactHistoryRef())
      {
         WALBERLA_CHECK_EQUAL(entry.second.getTangentialSpringDisplacement(), Vec3(real_c(entry.first * simulationSteps)), p);
      }
   }

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
