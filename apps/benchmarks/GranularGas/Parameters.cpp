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

void loadFromConfig(Parameters& params, const Config::BlockHandle& cfg)
{
   params.sorting = cfg.getParameter<std::string>("sorting", "none" );
   WALBERLA_LOG_INFO_ON_ROOT("sorting: " << params.sorting);
   
   params.spacing = cfg.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << params.spacing);
   
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
   
   params.path = cfg.getParameter<std::string>("path", "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("path: " << params.path);
   
   params.sqlFile = cfg.getParameter<std::string>("sqlFile", "benchmark.sqlite" );
   WALBERLA_LOG_INFO_ON_ROOT("sqlFile: " << params.sqlFile);
   
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
   
   stringProperties["path"] = params.path;
   
   stringProperties["sqlFile"] = params.sqlFile;
   
}

} //namespace walberla