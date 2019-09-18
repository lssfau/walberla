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
//! \file Parameters.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include "Parameters.h"

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

void loadFromConfig(Parameters& params, const Config::BlockHandle& cfg)
{
   params.sorting = cfg.getParameter<std::string>("sorting", "none" );
   WALBERLA_LOG_INFO_ON_ROOT("sorting: " << params.sorting);
   
   params.normal = cfg.getParameter<Vec3>("normal", Vec3(real_t(1.0), real_t(1.0), real_t(1.0)) );
   WALBERLA_LOG_INFO_ON_ROOT("normal: " << params.normal);
   
   params.spacing = cfg.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << params.spacing);
   
   params.shift = cfg.getParameter<Vec3>("shift", Vec3(real_t(0.1), real_t(0.1), real_t(0.1)) );
   WALBERLA_LOG_INFO_ON_ROOT("shift: " << params.shift);
   
   params.radius = cfg.getParameter<real_t>("radius", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << params.radius);
   
   params.bBarrier = cfg.getParameter<bool>("bBarrier", false );
   WALBERLA_LOG_INFO_ON_ROOT("bBarrier: " << params.bBarrier);
   
   params.storeNodeTimings = cfg.getParameter<bool>("storeNodeTimings", false );
   WALBERLA_LOG_INFO_ON_ROOT("storeNodeTimings: " << params.storeNodeTimings);
   
   params.checkSimulation = cfg.getParameter<bool>("checkSimulation", false );
   WALBERLA_LOG_INFO_ON_ROOT("checkSimulation: " << params.checkSimulation);
   
   params.numOuterIterations = cfg.getParameter<int64_t>("numOuterIterations", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("numOuterIterations: " << params.numOuterIterations);
   
   params.initialRefinementLevel = cfg.getParameter<int64_t>("initialRefinementLevel", 0 );
   WALBERLA_LOG_INFO_ON_ROOT("initialRefinementLevel: " << params.initialRefinementLevel);
   
   params.simulationSteps = cfg.getParameter<int64_t>("simulationSteps", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << params.simulationSteps);
   
   params.dt = cfg.getParameter<real_t>("dt", real_t(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << params.dt);
   
   params.visSpacing = cfg.getParameter<int64_t>("visSpacing", 1000 );
   WALBERLA_LOG_INFO_ON_ROOT("visSpacing: " << params.visSpacing);
   
   params.vtk_out = cfg.getParameter<std::string>("vtk_out", "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("vtk_out: " << params.vtk_out);
   
   params.sqlFile = cfg.getParameter<std::string>("sqlFile", "benchmark.sqlite" );
   WALBERLA_LOG_INFO_ON_ROOT("sqlFile: " << params.sqlFile);
   
   params.recalculateBlockLevelsInRefresh = cfg.getParameter<bool>("recalculateBlockLevelsInRefresh", false );
   WALBERLA_LOG_INFO_ON_ROOT("recalculateBlockLevelsInRefresh: " << params.recalculateBlockLevelsInRefresh);
   
   params.alwaysRebalanceInRefresh = cfg.getParameter<bool>("alwaysRebalanceInRefresh", true );
   WALBERLA_LOG_INFO_ON_ROOT("alwaysRebalanceInRefresh: " << params.alwaysRebalanceInRefresh);
   
   params.reevaluateMinTargetLevelsAfterForcedRefinement = cfg.getParameter<bool>("reevaluateMinTargetLevelsAfterForcedRefinement", false );
   WALBERLA_LOG_INFO_ON_ROOT("reevaluateMinTargetLevelsAfterForcedRefinement: " << params.reevaluateMinTargetLevelsAfterForcedRefinement);
   
   params.allowRefreshChangingDepth = cfg.getParameter<bool>("allowRefreshChangingDepth", false );
   WALBERLA_LOG_INFO_ON_ROOT("allowRefreshChangingDepth: " << params.allowRefreshChangingDepth);
   
   params.allowMultipleRefreshCycles = cfg.getParameter<bool>("allowMultipleRefreshCycles", false );
   WALBERLA_LOG_INFO_ON_ROOT("allowMultipleRefreshCycles: " << params.allowMultipleRefreshCycles);
   
   params.checkForEarlyOutInRefresh = cfg.getParameter<bool>("checkForEarlyOutInRefresh", true );
   WALBERLA_LOG_INFO_ON_ROOT("checkForEarlyOutInRefresh: " << params.checkForEarlyOutInRefresh);
   
   params.checkForLateOutInRefresh = cfg.getParameter<bool>("checkForLateOutInRefresh", true );
   WALBERLA_LOG_INFO_ON_ROOT("checkForLateOutInRefresh: " << params.checkForLateOutInRefresh);
   
   params.regridMin = cfg.getParameter<uint_t>("regridMin", uint_c(100) );
   WALBERLA_LOG_INFO_ON_ROOT("regridMin: " << params.regridMin);
   
   params.regridMax = cfg.getParameter<uint_t>("regridMax", uint_c(1000) );
   WALBERLA_LOG_INFO_ON_ROOT("regridMax: " << params.regridMax);
   
   params.maxBlocksPerProcess = cfg.getParameter<int>("maxBlocksPerProcess", int_c(1000) );
   WALBERLA_LOG_INFO_ON_ROOT("maxBlocksPerProcess: " << params.maxBlocksPerProcess);
   
   params.baseWeight = cfg.getParameter<real_t>("baseWeight", real_t(10.0) );
   WALBERLA_LOG_INFO_ON_ROOT("baseWeight: " << params.baseWeight);
   
   params.metisipc2redist = cfg.getParameter<real_t>("metisipc2redist", real_t(1000.0) );
   WALBERLA_LOG_INFO_ON_ROOT("metisipc2redist: " << params.metisipc2redist);
   
   params.LBAlgorithm = cfg.getParameter<std::string>("LBAlgorithm", "Hilbert" );
   WALBERLA_LOG_INFO_ON_ROOT("LBAlgorithm: " << params.LBAlgorithm);
   
   params.metisAlgorithm = cfg.getParameter<std::string>("metisAlgorithm", "PART_GEOM_KWAY" );
   WALBERLA_LOG_INFO_ON_ROOT("metisAlgorithm: " << params.metisAlgorithm);
   
   params.metisWeightsToUse = cfg.getParameter<std::string>("metisWeightsToUse", "BOTH_WEIGHTS" );
   WALBERLA_LOG_INFO_ON_ROOT("metisWeightsToUse: " << params.metisWeightsToUse);
   
   params.metisEdgeSource = cfg.getParameter<std::string>("metisEdgeSource", "EDGES_FROM_EDGE_WEIGHTS" );
   WALBERLA_LOG_INFO_ON_ROOT("metisEdgeSource: " << params.metisEdgeSource);
   
}

void saveToSQL(const Parameters& params,
               std::map< std::string, walberla::int64_t >& integerProperties,
               std::map< std::string, double >&            realProperties,
               std::map< std::string, std::string >&       stringProperties )
{
   stringProperties["sorting"] = params.sorting;
   
   
   realProperties["spacing"] = double_c(params.spacing);
   
   
   realProperties["radius"] = double_c(params.radius);
   
   
   
   
   integerProperties["numOuterIterations"] = params.numOuterIterations;
   
   integerProperties["initialRefinementLevel"] = params.initialRefinementLevel;
   
   integerProperties["simulationSteps"] = params.simulationSteps;
   
   realProperties["dt"] = double_c(params.dt);
   
   integerProperties["visSpacing"] = params.visSpacing;
   
   stringProperties["vtk_out"] = params.vtk_out;
   
   stringProperties["sqlFile"] = params.sqlFile;
   
   
   
   
   
   
   
   
   
   
   
   realProperties["baseWeight"] = double_c(params.baseWeight);
   
   realProperties["metisipc2redist"] = double_c(params.metisipc2redist);
   
   stringProperties["LBAlgorithm"] = params.LBAlgorithm;
   
   stringProperties["metisAlgorithm"] = params.metisAlgorithm;
   
   stringProperties["metisWeightsToUse"] = params.metisWeightsToUse;
   
   stringProperties["metisEdgeSource"] = params.metisEdgeSource;
   
}

} //namespace mesa_pd
} //namespace walberla