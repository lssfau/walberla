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
//! \file BodyStatistics.h
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <blockforest/StructuredBlockForest.h>

#include <core/DataTypes.h>
#include <core/math/DistributedSample.h>
#include <core/math/Sample.h>

#include <iostream>

namespace walberla {
namespace pe {

class BodyStatistics
{
public:

   BodyStatistics( const shared_ptr<BlockStorage>& blockStorage, const BlockDataID & bodyStorageID )
      : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID )
   { }

   void operator()();

   void toStream( std::ostream & os ) const;

   uint_t numBodies() const { return numBodies_; }
   uint_t numShadowCopies() const { return numShadowCopies_; }
   real_t minVelocity() const { return velocitySample_.min(); }
   real_t maxVelocity() const { return velocitySample_.max(); }
   real_t avgVelocity() const { return velocitySample_.avg(); }
   real_t totalMass() const { return massSample_.sum(); }

private:
   const shared_ptr<BlockStorage> blockStorage_;
   const BlockDataID              bodyStorageID_;

   math::Sample localBodiesBlockSample_, shadowBodiesBlockSample_;
   math::Sample localBodiesProcessSample_, shadowBodiesProcessSample_;
   uint_t numBodies_;
   uint_t numShadowCopies_;

   math::DistributedSample velocitySample_;
   math::DistributedSample massSample_;
};

inline std::ostream & operator<<( std::ostream & os, const BodyStatistics & bodyStatistics )
{
   bodyStatistics.toStream( os );
   return os;
}



} // namespace pe
} // namespace walberla

