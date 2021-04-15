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
//! \file BodyStatistics.cpp
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "BodyStatistics.h"

#include <core/logging/Logging.h>
#include <core/mpi/Reduce.h>

#include <pe/Types.h>
#include <pe/rigidbody/BodyStorage.h>
#include <pe/rigidbody/Sphere.h>

namespace walberla {
namespace pe {
   
void BodyStatistics::operator()()
{
   numBodies_ = 0;
   numShadowCopies_ = 0;
   localBodiesBlockSample_.clear();
   shadowBodiesBlockSample_.clear();
   localBodiesProcessSample_.clear();
   shadowBodiesProcessSample_.clear();
   velocitySample_.clear();
   massSample_.clear();

   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      auto storage = blockIt->getData< pe::Storage >( bodyStorageID_ );
      const pe::BodyStorage & localStorage  = (*storage)[0];
      const pe::BodyStorage & shadowStorage = (*storage)[1];

      numBodies_       += localStorage.size();
      numShadowCopies_ += shadowStorage.size();

      localBodiesBlockSample_.castToRealAndInsert( localStorage.size() );
      shadowBodiesBlockSample_.castToRealAndInsert( shadowStorage.size() );

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt )
      {
         velocitySample_.insert( bodyIt->getLinearVel().length() );
         massSample_.insert( bodyIt->getMass() );
      }
   }

   localBodiesProcessSample_.castToRealAndInsert( numBodies_ );
   shadowBodiesProcessSample_.castToRealAndInsert( numShadowCopies_ );

   localBodiesBlockSample_.mpiAllGather();
   shadowBodiesBlockSample_.mpiAllGather();
   localBodiesProcessSample_.mpiAllGather();
   shadowBodiesProcessSample_.mpiAllGather();
   mpi::allReduceInplace( numBodies_, mpi::SUM );
   mpi::allReduceInplace( numShadowCopies_, mpi::SUM );

   velocitySample_.mpiAllGather();
   massSample_.mpiAllGather();
}

void BodyStatistics::toStream( std::ostream & os ) const
{
   os << "Number of bodies:           " << numBodies_ << "\n" <<
         "Number of shadow copies:    " << numShadowCopies_ << "\n" <<
         "Bodies on blocks:           " << localBodiesBlockSample_.format() << "\n" <<
         "Shadow copies on blocks:    " << shadowBodiesBlockSample_.format() << "\n" <<
         "Bodies on processes:        " << localBodiesProcessSample_.format() << "\n" <<
         "Shadow copies on processes: " << shadowBodiesProcessSample_.format() << "\n" <<
         "Velocity:                   " << velocitySample_.format() << "\n" <<
         "Mass:                       " << massSample_.format();
}

} // namespace pe
} // namespace walberla

