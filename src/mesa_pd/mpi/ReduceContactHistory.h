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
//! \file ReduceContactHistory.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ContactHistory.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/kernel/ParticleSelector.h>

#include <mesa_pd/mpi/notifications/ContactHistoryNotification.h>
#include <mesa_pd/mpi/ReduceProperty.h>

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

/**
 * Reduces all contact history from ghost particles to the master particle.
 *
 * \attention
 * This kernel *only* reduces the contact history. You have to manually
 * broadcast it to ghost particles if needed. This can be done by either
 * using BroadcastProperty or SyncNextNeighbors.
 *
 * \par Usage:
 * The contact history aggregated by a call to this kernel is available in
 * oldContactHistory. Use this data to write the current contact history
 * to newContactHistory. The next call to this kernel will delete all
 * contact data which is not stored in newContactHistory!
 *
 * \pre
 * - up to date contact information has to be stored in newContactHistory
 * - only contact information of local contacts should be stored
 * - data in oldContactHistory will be overwritten
 *
 * \post
 * - oldContactHistory for the master particle contains the reduced data
 * - oldContactHistory for ghost particles is empty
 * - newContactHistory is empty for master as well as ghost particles
 *
 * \internal
 * Two contact histories are needed to decided which contacts are still persistent and
 * which have been separated. OldContactHistory stores all contact information from the last
 * time step. This can be used to populate newContactHistory. Processes should only populate
 * newContactHistory if they have a persisting contact! During the reduce operation newContactHistory
 * will be reduced effectively eliminating all contacts that have separated. Subsequently it
 * is swapped with oldContactHistory to create the initial state. This is oldContactHistory
 * has all information, newContactHistory is empty.
 *
 * \ingroup mesa_pd_mpi
 */
class ReduceContactHistory
{
public:
   void operator()(data::ParticleStorage& ps) const;
private:
   ReduceProperty RP_;

   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();
};

void ReduceContactHistory::operator()(data::ParticleStorage& ps) const
{
   //no need to reduce if run with only one process
   if (numProcesses_ != 1)
   {
      RP_.operator()<ContactHistoryNotification>(ps);
   }

   const auto size = ps.size();
   for (size_t idx = 0; idx < size; ++idx)
   {
      if (!data::particle_flags::isSet( ps.getFlags(idx), data::particle_flags::GHOST) )
      {
         std::swap(ps.getOldContactHistoryRef(idx), ps.getNewContactHistoryRef(idx));
         ps.getNewContactHistoryRef(idx).clear();
      } else
      {
         ps.getOldContactHistoryRef(idx).clear();
         ps.getNewContactHistoryRef(idx).clear();
      }
   }
}

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla