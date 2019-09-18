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
//! \file   DetectAndStoreCollisions.cpp
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================


/** Test Collision Detection and Insertion of contacts into the contact storage */

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ContactStorage.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/kernel/DetectAndStoreContacts.h>
#include <mesa_pd/domain/InfiniteDomain.h>
#include <mesa_pd/kernel/ParticleSelector.h>

#include <mesa_pd/data/ParticleAccessor.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

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

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);

   auto  smallSphere = ss->create<data::Sphere>( real_t(1.2) );

   ss->shapes[smallSphere]->updateMassAndInertia(real_t(1));

   domain::InfiniteDomain domain;

   // Create four slightly overlapping spheres in a row (located at x=0,2,4,6)
   for (int i = 0; i < 8; i+=2)
   {
      auto p                       = ps->create();
      p->getPositionRef()          = Vec3(real_t(i), real_t(0), real_t(0));
      p->getShapeIDRef()           = smallSphere;
      p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
      p->getTypeRef()              = 0;
   }


   Vec3 normal(-1,0,0);
   auto dist= real_t(-0.4);

   // Create Contact Storage cs
   data::ContactStorage cs(100);
   cs.clear();

   // Perform Collision detection (call kernel, that stores contacts into cs)
   kernel::DetectAndStoreContacts detectAndStore(cs);
   ps->forEachParticlePairHalf(false, kernel::ExcludeInfiniteInfinite(), accessor, detectAndStore, accessor, domain);

   // Check if all three intersections were found
   size_t contactCount = cs.size();
   WALBERLA_CHECK_EQUAL(contactCount, 3);

   // Check the contacts with the for each Contact loop.
   // Set back contact count to 0 to now count the loop iterations.
   contactCount = 0;

   cs.forEachContact(false, kernel::SelectAll(), cs, [&normal, &dist, &contactCount](size_t idx, data::ContactStorage &css){
      WALBERLA_CHECK_FLOAT_EQUAL(css.getNormal(idx), normal);
      WALBERLA_CHECK_FLOAT_EQUAL(css.getPosition(idx), Vec3(real_t(2*idx+1), real_t(0), real_t(0)));
      WALBERLA_CHECK_FLOAT_EQUAL(css.getDistance(idx), dist);
      contactCount++;
   }
   ,cs);

   WALBERLA_CHECK_EQUAL(contactCount, 3);

   WALBERLA_LOG_INFO("Insertion test with ContactStorage successful.");
   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
