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
//! \file InspectionProbe.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "mesa_pd/data/Flags.h"

namespace walberla {
namespace lbm_mesapd_coupling {


/*
 * Functionality that allows to monitor the state of a specific cell in the fluid-particle simulation.
 * Options to output the state (fluid with density and velocity, or particle) on screen or to file are available.
 * Also, the state of surrounding cells can be printed to screen.
 * Can generally be used to track down the source of instabilities in a specific cell which often has its origin locally in the direct neighborhood of the cell.
 *
 */
template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T >
class InspectionProbe
{

public:
   InspectionProbe(Vector3<real_t> probeLocation,
                   const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & pdfFieldID, const BlockDataID & boundaryHandlingID,
                   const BlockDataID & particleFieldID, const shared_ptr< ParticleAccessor_T > & ac,
                   bool printToScreen, bool printSurroundingState, std::string outputFileName)
   : probeLocation_(probeLocation), blocks_(blocks), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID), particleFieldID_(particleFieldID),
     ac_(ac), printToScreen_(printToScreen), printSurroundingState_(printSurroundingState), outputFileName_(outputFileName)
   {}

   void operator()(real_t & rho, Vector3<real_t> & velocity)
   {

      Cell probeCellGlobal = blocks_->getCell(probeLocation_);
      auto block = blocks_->getBlock(probeLocation_);
      if(block != nullptr)
      {
         auto pdfField = block->template getData<PdfField_T>(pdfFieldID_);
         auto boundaryHandling = block->template getData<BoundaryHandling_T>(boundaryHandlingID_);

         Cell probeCell;
         blocks_->transformGlobalToBlockLocalCell( probeCell, *block, probeCellGlobal );

         if(printSurroundingState_) printStatusOfCellSurrounding(probeCell, *block);

         rho = pdfField->getDensityAndVelocity(velocity,probeCell);

         bool isFluid = boundaryHandling->isDomain(probeCell);

         printToScreen(isFluid, rho, velocity);
         writeToFile(isFluid, rho, velocity);

      }
   }

   void setPosition(Vector3<real_t> probeLocation)
   {
      probeLocation_ = probeLocation;
   }

private:

   void printStatusOfCellSurrounding(Cell centerCell, const IBlock& block)
   {
      auto pdfField = block.getData<PdfField_T>(pdfFieldID_);
      auto boundaryHandling = block.getData<BoundaryHandling_T>(boundaryHandlingID_);
      auto particleField = block.getData<lbm_mesapd_coupling::ParticleField_T>(particleFieldID_);

      std::stringstream outputString;

      for(cell_idx_t z = cell_idx_t(-1); z <= cell_idx_t(1); ++z)
      {
         for(cell_idx_t y = cell_idx_t(1); y >= cell_idx_t(-1); --y)
         {
            outputString << "*------------------------------------------------------------------------------------------------------------------------------------------------------------------------*\n";
            for(cell_idx_t x = cell_idx_t(-1); x <= cell_idx_t(1); ++x)
            {
               auto cell = centerCell + Cell(x,y,z);
               auto cellCenter = blocks_->getBlockLocalCellCenter(block, cell);
               real_t rho;
               Vector3<real_t> velocity;
               bool isFluid=false;
               bool isFixed=false;
               bool isGlobal=false;
               if(boundaryHandling->isDomain(cell))
               {
                  isFluid = true;
                  rho = pdfField->getDensityAndVelocity(velocity,cell);
               } else{
                  isFluid = false;
                  rho = real_t(0);
                  auto particleIdx = ac_->uidToIdx(particleField->get( cell ));
                  velocity = mesa_pd::getVelocityAtWFPoint(particleIdx, *ac_, cellCenter );
                  if(isSet(ac_->getFlags(particleIdx), mesa_pd::data::particle_flags::FIXED)) isFixed = true;
                  if(isSet(ac_->getFlags(particleIdx), mesa_pd::data::particle_flags::GLOBAL)) isGlobal = true;
               }
               outputString << std::setprecision(5) << "| " << cellCenter << " (" << (isFluid?"F":(std::string("P")+(isFixed?"+F":"")+(isGlobal?"+G":""))) << ") " << rho << " " << velocity << " ";
            }
            outputString << "|\n";
         }
         outputString << "*------------------------------------------------------------------------------------------------------------------------------------------------------------------------*\n\n";
      }
      WALBERLA_LOG_INFO(outputString.str());
   }

   void printToScreen(bool isFluid, real_t rho, Vector3<real_t> velocity)
   {
      if(printToScreen_) WALBERLA_LOG_INFO("Values in probe at position " << probeLocation_ << " on rank " << mpi::MPIManager::instance()->rank() << ": " << (isFluid?"F":"P") << ", rho = " << rho << ", velocity = " << velocity);
   }

   void writeToFile(bool isFluid, real_t rho, Vector3<real_t> velocity)
   {
      if(!outputFileName_.empty())
      {
         std::ofstream file;
         file.open( outputFileName_.c_str(), std::ofstream::app);

         file << isFluid << " " << rho << " " << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
         file.close();
      }
   }

   Vector3<real_t> probeLocation_;
   shared_ptr< StructuredBlockStorage > blocks_;
   const ConstBlockDataID pdfFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
   const ConstBlockDataID particleFieldID_;
   shared_ptr< ParticleAccessor_T > ac_;
   bool printToScreen_;
   bool printSurroundingState_;
   std::string outputFileName_;

};

} // namespace lbm_mesapd_coupling
} // namespace walberla
