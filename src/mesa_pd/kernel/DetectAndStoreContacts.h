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
//! \file DetectAndStoreContacts.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/mpi/ContactFilter.h>

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ContactStorage.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/data/Flags.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Kernel which performes collision detection on a pair of two particles
 * and inserts the contact (if existent) into the contact storage passed in the constructor.
 * Call this kernel on each particle pair to perform contact detection and insert each contact in the contact
 * storage.
 * \ingroup mesa_pd_kernel
 */
class DetectAndStoreContacts
{

   public:
   explicit DetectAndStoreContacts(data::ContactStorage &cs) : cs_(cs) {}

   /**
    * Compute if particle 1 and particle 2 collide. If so, insert the contact into the contact storage.
    * \param idx1 The index of particle 1
    * \param idx2 The index of particle 2
    * \param ac The accessor used to access the values of the particles.
    * \param domain The domain of the block (used to decide if the contact has to be treated by this process)
    * \param acd The collision detection to be used. Default parameter: AnalyticContactDetection()
    */
   template <typename Accessor>
   void operator()(size_t idx1, size_t idx2, Accessor &ac, const domain::IDomain& domain, collision_detection::AnalyticContactDetection acd = collision_detection::AnalyticContactDetection());
   private:
   data::ContactStorage& cs_;
};



template <typename Accessor>
inline void DetectAndStoreContacts::operator()(size_t idx1, size_t idx2, Accessor &ac, const domain::IDomain& domain, collision_detection::AnalyticContactDetection acd)
{
   using namespace data::particle_flags;
   kernel::DoubleCast double_cast;
   mpi::ContactFilter contact_filter;
   if (double_cast(idx1, idx2, ac, acd, ac ))
   {
      if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
      {
         auto c = cs_.create();
         c->setId1(acd.getIdx1());
         c->setId2(acd.getIdx2());
         c->setDistance(acd.getPenetrationDepth());
         c->setNormal(acd.getContactNormal());
         c->setPosition(acd.getContactPoint());
      }
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla