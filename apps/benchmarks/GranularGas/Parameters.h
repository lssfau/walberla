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
//! \file Parameters.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <core/config/Config.h>
#include <mesa_pd/data/DataTypes.h>

#include <string>

namespace walberla {
namespace mesa_pd {

struct Parameters
{
   std::string sorting = "none";
   Vec3 normal = Vec3(real_t(1.0), real_t(1.0), real_t(1.0));
   real_t spacing = real_t(1.0);
   Vec3 shift = Vec3(real_t(0.1), real_t(0.1), real_t(0.1));
   real_t radius = real_t(0.5);
   bool bBarrier = false;
   bool storeNodeTimings = false;
   bool checkSimulation = false;
   int64_t numOuterIterations = 10;
   int64_t initialRefinementLevel = 0;
   int64_t simulationSteps = 10;
   real_t dt = real_t(0.01);
   int64_t visSpacing = 1000;
   std::string vtk_out = "vtk_out";
   std::string sqlFile = "benchmark.sqlite";
   bool recalculateBlockLevelsInRefresh = false;
   bool alwaysRebalanceInRefresh = true;
   bool reevaluateMinTargetLevelsAfterForcedRefinement = false;
   bool allowRefreshChangingDepth = false;
   bool allowMultipleRefreshCycles = false;
   bool checkForEarlyOutInRefresh = true;
   bool checkForLateOutInRefresh = true;
   uint_t regridMin = uint_c(100);
   uint_t regridMax = uint_c(1000);
   int maxBlocksPerProcess = int_c(1000);
   real_t baseWeight = real_t(10.0);
   real_t metisipc2redist = real_t(1000.0);
   std::string LBAlgorithm = "Hilbert";
   std::string metisAlgorithm = "PART_GEOM_KWAY";
   std::string metisWeightsToUse = "BOTH_WEIGHTS";
   std::string metisEdgeSource = "EDGES_FROM_EDGE_WEIGHTS";
};

void loadFromConfig(Parameters& params,
                    const Config::BlockHandle& cfg);

void saveToSQL(const Parameters& params,
               std::map< std::string, walberla::int64_t >& integerProperties,
               std::map< std::string, double >&            realProperties,
               std::map< std::string, std::string >&       stringProperties );

} //namespace mesa_pd
} //namespace walberla