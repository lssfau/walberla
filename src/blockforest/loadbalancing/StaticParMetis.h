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
//! \file StaticParMetis.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/SetupBlockForest.h"

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIWrapper.h"

#include <map>
#include <vector>

namespace walberla {
namespace blockforest {


class StaticLevelwiseParMetis
{
public:
   enum Algorithm    { PARMETIS_PART_GEOM_KWAY, PARMETIS_PART_KWAY };
   enum WeightsToUse { PARMETIS_NO_WEIGHTS = 0, PARMETIS_EDGE_WEIGHTS = 1, PARMETIS_VERTEX_WEIGHTS = 2, PARMETIS_BOTH_WEIGHTS = 3 };

   using BlockPair = std::pair<const SetupBlock *, const SetupBlock *>;
   using CommWeightFunction = std::function<void (const std::vector<BlockPair> &, std::vector<int64_t> &)>;

   StaticLevelwiseParMetis( const Algorithm algorithm = PARMETIS_PART_GEOM_KWAY )
      : algorithm_( algorithm ), weightsToUse_( PARMETIS_VERTEX_WEIGHTS ) {}

   StaticLevelwiseParMetis( const CommWeightFunction & commWeightFunction, 
                            const Algorithm algorithm = PARMETIS_PART_GEOM_KWAY,
                            const WeightsToUse weightsToUse = PARMETIS_BOTH_WEIGHTS )
      : algorithm_( algorithm ), weightsToUse_( weightsToUse ), commWeightFunction_( commWeightFunction ) {}

   uint_t operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t perProcessMemoryLimit ) const;

   bool edgeWeightsUsed()   const { return ( weightsToUse_ == PARMETIS_EDGE_WEIGHTS   ) || ( weightsToUse_ == PARMETIS_BOTH_WEIGHTS ); }
   bool vertexWeightsUsed() const { return ( weightsToUse_ == PARMETIS_VERTEX_WEIGHTS ) || ( weightsToUse_ == PARMETIS_BOTH_WEIGHTS ); }

   static Algorithm    stringToAlgorithm( std::string s );
   static WeightsToUse stringToWeightsToUse( std::string s );

protected:
   Algorithm algorithm_;
   WeightsToUse weightsToUse_;
   CommWeightFunction commWeightFunction_;
};



} // namespace blockforest
} // namespace walberla
