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

#pragma once

#include "core/config/Config.h"

#include <string>

namespace walberla {
namespace mesa_pd {

struct FilmSpecimen
{
   real_t sizeX = 0; //Specimen length (x direction)
   real_t sizeY = 0; //Specimen length (y direction)
   real_t sizeZ = 0; //Specimen length (z direction)
   bool oopp = false; //Specimen out-of-plane periodicity (0 - film, 1 - cnt material)
   uint_t numBlocksX = 1; //Number of blocks in x direction
   uint_t numBlocksY = 1; //Number of blocks in y direction
   uint_t numBlocksZ = 1; //Number of blocks in z direction
   real_t min_OOP = 0; //Out-of-plane angle minimum
   real_t max_OOP = 0; //Out-of-plane angle maximum
   int numCNTs = 0; //Number of CNTs
   int numSegs = 0; //Number of segments in a CNT
   real_t spacing = 0; //Segment half-spacing
   real_t localDamping = 0; //Local damping coefficient
   real_t viscousDamping = 0; //Viscous damping coefficient
   uint_t seed = 0; //random generator seed
   int vdW = 0; //type of vdW interaction model
   int simulationSteps = 0; //Relaxation duration
   int saveVTKEveryNthStep = 0; //timesteps between saving VTK outputs
   int saveEnergyEveryNthStep = 0; //timesteps between saving energies
   int saveConfEveryNthStep = 0; //timesteps between saving confs
   std::string vtkFolder = "."; //Folder for VTK files
   std::string energyFolder = "."; //Folder for energy files
   std::string confFolder = "."; //Folder for conf files
   bool useMPIIO = false; //Write a single file instead of one file per process
   std::string sqlFile = "cnt.sqlite"; //database file
   std::string initialConfigurationFile = ""; //restart from checkpoint
};

void loadFromConfig(FilmSpecimen& params,
                    const Config::BlockHandle& cfg);

} //namespace mesa_pd
} //namespace walberla