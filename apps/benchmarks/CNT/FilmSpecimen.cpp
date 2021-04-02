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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include "FilmSpecimen.h"
#include "TerminalColors.h"

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

void loadFromConfig(FilmSpecimen& params, const Config::BlockHandle& cfg)
{
   WALBERLA_LOG_INFO_ON_ROOT(GREEN << "Loading CNT film setup: " << RESET)
   params.sizeX = cfg.getParameter<real_t>("sizeX", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "sizeX: " << GREEN << params.sizeX << RESET);
   
   params.sizeY = cfg.getParameter<real_t>("sizeY", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "sizeY: " << GREEN << params.sizeY << RESET);
   
   params.sizeZ = cfg.getParameter<real_t>("sizeZ", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "sizeZ: " << GREEN << params.sizeZ << RESET);
   
   params.oopp = cfg.getParameter<bool>("oopp", false );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "oopp: " << GREEN << params.oopp << RESET);
   
   params.numBlocksX = cfg.getParameter<uint_t>("numBlocksX", 1 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "numBlocksX: " << GREEN << params.numBlocksX << RESET);
   
   params.numBlocksY = cfg.getParameter<uint_t>("numBlocksY", 1 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "numBlocksY: " << GREEN << params.numBlocksY << RESET);
   
   params.numBlocksZ = cfg.getParameter<uint_t>("numBlocksZ", 1 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "numBlocksZ: " << GREEN << params.numBlocksZ << RESET);
   
   params.min_OOP = cfg.getParameter<real_t>("min_OOP", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "min_OOP: " << GREEN << params.min_OOP << RESET);
   
   params.max_OOP = cfg.getParameter<real_t>("max_OOP", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "max_OOP: " << GREEN << params.max_OOP << RESET);
   
   params.numCNTs = cfg.getParameter<int>("numCNTs", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "numCNTs: " << GREEN << params.numCNTs << RESET);
   
   params.numSegs = cfg.getParameter<int>("numSegs", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "numSegs: " << GREEN << params.numSegs << RESET);
   
   params.spacing = cfg.getParameter<real_t>("spacing", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "spacing: " << GREEN << params.spacing << RESET);
   
   params.localDamping = cfg.getParameter<real_t>("localDamping", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "localDamping: " << GREEN << params.localDamping << RESET);
   
   params.viscousDamping = cfg.getParameter<real_t>("viscousDamping", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "viscousDamping: " << GREEN << params.viscousDamping << RESET);
   
   params.seed = cfg.getParameter<uint_t>("seed", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "seed: " << GREEN << params.seed << RESET);
   
   params.vdW = cfg.getParameter<int>("vdW", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "vdW: " << GREEN << params.vdW << RESET);
   
   params.simulationSteps = cfg.getParameter<int>("simulationSteps", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "simulationSteps: " << GREEN << params.simulationSteps << RESET);
   
   params.saveVTKEveryNthStep = cfg.getParameter<int>("saveVTKEveryNthStep", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "saveVTKEveryNthStep: " << GREEN << params.saveVTKEveryNthStep << RESET);
   
   params.saveEnergyEveryNthStep = cfg.getParameter<int>("saveEnergyEveryNthStep", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "saveEnergyEveryNthStep: " << GREEN << params.saveEnergyEveryNthStep << RESET);
   
   params.saveConfEveryNthStep = cfg.getParameter<int>("saveConfEveryNthStep", 0 );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "saveConfEveryNthStep: " << GREEN << params.saveConfEveryNthStep << RESET);
   
   params.vtkFolder = cfg.getParameter<std::string>("vtkFolder", "." );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "vtkFolder: " << GREEN << params.vtkFolder << RESET);
   
   params.energyFolder = cfg.getParameter<std::string>("energyFolder", "." );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "energyFolder: " << GREEN << params.energyFolder << RESET);
   
   params.confFolder = cfg.getParameter<std::string>("confFolder", "." );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "confFolder: " << GREEN << params.confFolder << RESET);
   
   params.useMPIIO = cfg.getParameter<bool>("useMPIIO", false );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "useMPIIO: " << GREEN << params.useMPIIO << RESET);
   
   params.sqlFile = cfg.getParameter<std::string>("sqlFile", "cnt.sqlite" );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "sqlFile: " << GREEN << params.sqlFile << RESET);
   
   params.initialConfigurationFile = cfg.getParameter<std::string>("initialConfigurationFile", "" );
   WALBERLA_LOG_INFO_ON_ROOT(YELLOW << "initialConfigurationFile: " << GREEN << params.initialConfigurationFile << RESET);
   
}

} //namespace mesa_pd
} //namespace walberla

