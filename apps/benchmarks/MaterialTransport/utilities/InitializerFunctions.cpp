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
//! \file InitializerFunctions.h
//! \author Ravi Ayyala Somayajula <ravi.k.ayyala@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/all.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

// #include "python_coupling/DictWrapper.h"
#include "GeneralInfoHeader.h"
#include "InitializerFunctions.h"
#pragma once

namespace walberla
{

void initConcentrationField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& ConcentrationFieldID,
                            const math::AABB& domainAABB, Vector3< uint_t > domainSize)
{
   const real_t radius = real_c(domainSize[1] / 10);
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      WALBERLA_LOG_INFO_ON_ROOT("failed here")
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(ConcentrationField, {
         Cell globalCell;
         const auto cellAABB = blocks->getBlockLocalCellAABB(block, Cell(x, y, z));
         auto cellCenter     = cellAABB.center();
         const real_t posX   = x; // cellCenter[0];
         const real_t posY   = y; // cellCenter[1];
         const real_t posZ   = z; // cellCenter[2];

         real_t distance =
            real_c(sqrt(pow((domainAABB.center()[0] - posX), 2) + pow((domainAABB.center()[1] - posY), 2) +
                        pow((domainAABB.center()[2] - posZ), 2)));

         /*if (distance <= radius) { ConcentrationField->get(x, y, z) = real_t(0.0); }
         else { ConcentrationField->get(x, y, z) = real_t(0.0); }*/
         ConcentrationField->get(x, y, z) = real_t(0.0);
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

} // initConcentrationField

void initConcentrationFieldGaussian(const shared_ptr< StructuredBlockStorage >& blocks,
                                    BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                    Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                    const Vector3< real_t > uInflow, const Vector3< real_t > x_0)
{
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      Block& b                = dynamic_cast< Block& >(block);
      uint_t level            = b.getLevel();
      CellInterval xyz        = ConcentrationField->xyzSize();
      double sum=0;
      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         // WALBERLA_LOG_INFO("posx " << pos[0]);
         ConcentrationField->get(*cellIt) = std::exp(
            -(std::pow((pos[0] - x_0[0]), 2) + std::pow((pos[1] - x_0[1]), 2) + std::pow((pos[2]-x_0[2]),2)) /
            (2 * sigma_0 * sigma_0));
         //ConcentrationField->get(*cellIt) = std::max(ConcentrationField->get(*cellIt), 1e-15);

      }
   }
}

void initConcentrationFieldSinusoidal(const shared_ptr< StructuredBlockStorage >& blocks,
                                    BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                    Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                    const Vector3< real_t > uInflow, const Vector3< real_t > x_0,const real_t dx,const real_t dt)
{
   const real_t pi = M_PI;
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      Block& b                = dynamic_cast< Block& >(block);
      uint_t level            = b.getLevel();
      CellInterval xyz        = ConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         ConcentrationField->get(*cellIt) = std::sin(pi*pos[0]*dx)*std::sin(pi*pos[1]*dx)*std::sin(pi*pos[2]*dx) +1;
         ConcentrationField->get(*cellIt) = std::max(ConcentrationField->get(*cellIt), 1e-15);
      }
   }
}

void initConcentrationFieldPacket(const shared_ptr< StructuredBlockStorage >& blocks,
                                      BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                      Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                      const Vector3< real_t > uInflow, const Vector3< real_t > x_0,const real_t dx,const real_t dt,const real_t diffusivity)
{
   const real_t pi = M_PI;
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      Block& b                = dynamic_cast< Block& >(block);
      uint_t level            = b.getLevel();
      CellInterval xyz        = ConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         if(pos[0]*dx < x_0[0]*dx + 0.5*dx && pos[0]*dx > x_0[0]*dx - 0.5*dx && pos[1]*dx < x_0[1]*dx + 0.5*dx && pos[1]*dx > x_0[1]*dx - 0.5*dx && pos[2]*dx < x_0[2]*dx + 0.5*dx && pos[2]*dx > x_0[2]*dx - 0.5*dx){
            ConcentrationField->get(*cellIt) = 1/std::sqrt(4*pi*diffusivity*dt) + 1;
            WALBERLA_LOG_INFO("True ");
         }
         else{
            ConcentrationField->get(*cellIt) = real_c(1);
         }
      }
   }
}

void initFluidField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,
                    const Vector3< real_t > uInflow,Vector3< uint_t > domainSize)
{
   for (auto& block : *blocks)
   {

      auto FluidVelocityField = block.getData< VelocityField_fluid_T >(FluidFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(FluidVelocityField, {
         Cell globalCell;
         const auto cellAABB                 = blocks->getBlockLocalCellAABB(block, Cell(x, y, z));
         FluidVelocityField->get(x, y, z, 0) = 0;//uInflow[0];
         FluidVelocityField->get(x, y, z, 1) = 0;//uInflow[1];
         if(domainSize[2] != 1)
         {
            FluidVelocityField->get(x, y, z, 2) = 0;
         }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }
}

void initFluidFieldPoiseuille(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,
                    const Vector3< real_t > uInflow,Vector3< uint_t > domainSize)
{
   for (auto& block : *blocks)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Velocity init reached here");
      //WALBERLA_LOG_INFO_ON_ROOT("y height is " << domainSize[1]);
      auto FluidVelocityField = block.getData< VelocityField_fluid_T >(FluidFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(FluidVelocityField, {
         Cell globalCell;
         const auto cellAABB                 = blocks->getBlockLocalCellAABB(block, Cell(x, y, z));
         FluidVelocityField->get(x, y, z, 0) = uInflow[0]*(1 - std::pow((y/real_c(domainSize[1])),2));
         FluidVelocityField->get(x, y, z, 1) = 0;

      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }
}



void analyticalSolGaussian(const shared_ptr< StructuredBlockStorage >& blocks,
                           BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                           Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                           const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time)
{
   real_t sigma_D2 = (2 * diffusivity * time);

   for (auto& block : *blocks)
   {
      auto AnalyticalConcentrationField = block.getData< DensityField_concentration_T >(AnalyticalConcentrationFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      uint_t level                      = b.getLevel();
      CellInterval xyz                  = AnalyticalConcentrationField->xyzSize();
      double sum = 0.0;
      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         // WALBERLA_LOG_INFO("posx " << pos[0]);
         real_t prefactor = (sigma_0 * sigma_0) / (sigma_0 * sigma_0 + sigma_D2);
         //real_t prefactor = (std::pow(sigma_0,3)) /(std::pow((sigma_0 * sigma_0 + sigma_D2),1.5));

         double x_eff = pos[0] - x_0[0] - uInflow[0] * time;
         double y_eff = pos[1] - x_0[1] - uInflow[1] * time;
         /*if(x_eff < -advection_period*domainSize[0]/2){
            x_eff = x_eff + advection_period*domainSize[0];
         }
         if(y_eff < -advection_period*domainSize[1]/2){
            y_eff = y_eff + advection_period*domainSize[1];
         }
         if(x_eff > advection_period*domainSize[0]/2){
            x_eff = x_eff - advection_period*domainSize[0];
         }
         if(y_eff > advection_period*domainSize[1]/2){
            y_eff = y_eff - advection_period*domainSize[1];
         }*/


         /*AnalyticalConcentrationField->get(*cellIt) =
            prefactor * std::exp(-(std::pow((x_eff + advection_period*domainSize[0]), 2) +
                                   std::pow((y_eff + advection_period*domainSize[1]),
                                            2)) /
                                 (2 * (sigma_0 * sigma_0 + sigma_D2)));*/

      }
   } /*+ std::pow((pos[2] - x_0[2] - uInflow[2]*time),2)*/
}

void analyticalSolSinusoidal(const shared_ptr< StructuredBlockStorage >& blocks,
                           BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                           Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                           const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time,const real_t dx,const real_t dt)
{
   const real_t pi = M_PI;
   real_t sigma_D2 = (2 * diffusivity * time);

   for (auto& block : *blocks)
   {
      auto AnalyticalConcentrationField = block.getData< DensityField_concentration_T >(AnalyticalConcentrationFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      uint_t level                      = b.getLevel();
      CellInterval xyz                  = AnalyticalConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         AnalyticalConcentrationField->get(*cellIt) = (std::sin(pi*(pos[0]*dx - 2.5*time*dt))*std::sin(pi*(pos[1]*dx - 2.5*time*dt))*std::sin(pi*(pos[2]*dx - 2.5*time*dt)))*std::exp(-3*pi*pi*time*dt*diffusivity) + 1;

      }
   }
}

void analyticalSolPacket(const shared_ptr< StructuredBlockStorage >& blocks,
                             BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                             Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                             const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time,const real_t dx,const real_t dt)
{
   const real_t pi = M_PI;
   real_t sigma_D2 = (2 * diffusivity * time);

   for (auto& block : *blocks)
   {
      auto AnalyticalConcentrationField = block.getData< DensityField_concentration_T >(AnalyticalConcentrationFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      uint_t level                      = b.getLevel();
      CellInterval xyz                  = AnalyticalConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);
         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         double sum = 0.0;
         if(time*dt == 0){
            const real_t expression = (pos[0] - x_0[0] - std::floor((pos[0] - x_0[0])/domainSize[0])*domainSize[0]) + (pos[1] - x_0[1] - std::floor((pos[1] - x_0[1])/domainSize[1])*domainSize[1]) + (pos[2] - x_0[2] - std::floor((pos[2] - x_0[2])/domainSize[2])*domainSize[2]);
            if(expression <= abs(1e-6)){AnalyticalConcentrationField->get(*cellIt) = real_c(1);}
            else{AnalyticalConcentrationField->get(*cellIt) = real_c(0);}
         }
         /*else{
            const real_t k = std::floor((pos[0] - x_0[0])/domainSize[0]);
            AnalyticalConcentrationField->get(*cellIt) = (1/(std::sqrt(4*pi*diffusivity*time*dt)))*std::exp(-((pos[0]*dx - x_0[0]*dx - 2.5*time*dt + k*domainSize[0]*dx) + (pos[1]*dx - x_0[1]*dx - 2.5*time*dt + k*domainSize[1]*dx) +(pos[2]*dx - x_0[2]*dx - 2.5*time*dt + k*domainSize[2]*dx))/(4*diffusivity*time*dt))+ 1;
         }*/

         // Sum over periodic images
         else
         {


               sum += (1/(std::sqrt(4*pi*diffusivity*time*dt)))*std::exp(-((pos[0]*dx - x_0[0]*dx - 2.5*time*dt + std::floor((pos[0] - x_0[0])/domainSize[0])*domainSize[0]*dx) + (pos[1]*dx - x_0[1]*dx - 2.5*time*dt + std::floor((pos[0] - x_0[0])/domainSize[0])*domainSize[1]*dx) +(pos[2]*dx - x_0[2]*dx - 2.5*time*dt + std::floor((pos[0] - x_0[0])/domainSize[0])*domainSize[2]*dx))/(4*diffusivity*time*dt))+ 1;


            AnalyticalConcentrationField->get(*cellIt) = sum;
         }

      }
   }
}

std::vector< real_t > computeErrorL2(const shared_ptr< StructuredBlockStorage >& blocks,
                                     BlockDataID& NumericalSolFieldID, BlockDataID& AnalyticalSolFieldID,
                                     BlockDataID& ErrorFieldID, const math::AABB& domainAABB)
{
   real_t Linf{ 0.0 };
   real_t L1{ 0.0 };
   real_t L2{ 0.0 };
   real_t analytical_squared{ 0.0 };
   uint_t cells{ 0 };
   std::vector< real_t > Errors{ 0, 0, 0 };

   for (auto block = blocks->begin(); block != blocks->end(); ++block)
   {
      auto numericalSolutionField  = block->getData< DensityField_concentration_T >(NumericalSolFieldID);
      auto analyticalSolutionField = block->getData< DensityField_concentration_T >(AnalyticalSolFieldID);
      auto ErrorField              = block->getData< DensityField_concentration_T >(ErrorFieldID);
      Block& b                     = dynamic_cast< Block& >(*block);
      uint_t level                 = b.getLevel();
      CellInterval xyz             = analyticalSolutionField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         real_t currErr = (numericalSolutionField->get(*cellIt) - analyticalSolutionField->get(*cellIt));

         // potentialErrorField ->get(x, y, z) = real_c(fabs(currErr));
         ErrorField->get(*cellIt) = fabs(currErr);
         L1 += std::abs(currErr);
         L2 += currErr * currErr;
         analytical_squared += analyticalSolutionField->get(*cellIt) * analyticalSolutionField->get(*cellIt);
         Linf = std::max(Linf, std::abs(currErr));
         cells += 1;
      }
   }
   mpi::allReduceInplace(L1, mpi::SUM);
   mpi::allReduceInplace(L2, mpi::SUM);
   mpi::allReduceInplace(analytical_squared, mpi::SUM);
   mpi::allReduceInplace(Linf, mpi::MAX);
   mpi::allReduceInplace(cells, mpi::SUM);
   Errors[0] = Linf;
   Errors[1] = L1 / cells;
   Errors[2] = std::sqrt(L2 / analytical_squared);
   return Errors;
}

std::vector< real_t > NusseltNumbers(const shared_ptr< StructuredBlockStorage >& blocks,
                             BlockDataID& ConcentrationFieldID,BlockDataID& VelocityFieldID, const math::AABB& domainAABB,
                             Vector3< uint_t > domainSize,const real_t delta_theta,const real_t ratio)
{
   std::vector< real_t > nusseltNumbers{ 0, 0, 0 };
   real_t nusseltx_domain = 0;
   real_t nusseltx_zero = 0; //Simpsons (real_c(128),real_c(0),ConcentrationFieldID,blocks);;
   real_t nusseltx_half = 0;
   const real_t Nx = real_c(domainSize[0]);
   const real_t Ny = real_c(domainSize[1]);

   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      auto VelocityField = block.getData< VelocityField_fluid_T >(VelocityFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      uint_t level                      = b.getLevel();
      CellInterval xyz                  = ConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);
         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         const cell_idx_t i = cellIt->x();
         const cell_idx_t j = cellIt->y();
         const cell_idx_t k = cellIt->z();
         //WALBERLA_LOG_INFO_ON_ROOT("i " << i << " " << "j " << j << " " << "k " << k);
         //const non_dim_temp = ConcentrationField->get(1,j,k) - 0.5;
         //const real_t differential_x_0 = 2*Nx*(ConcentrationField->get(0,j,k) - 1.0);
         //const real_t differential_x_half = Nx*(ConcentrationField->get(cell_idx_c(Nx/2),j,k) - ConcentrationField->get(cell_idx_c(Nx/2)-1,j,k));
         //const real_t differential_x_domain = Nx*(ConcentrationField->get(i+1,j,k) - ConcentrationField->get(i-1,j,k))/2;
         //const real_t differential_x_nondomain = Nx*(ConcentrationField->get(i+1,j,k) + ConcentrationField->get(i,j,k) - 2*ConcentrationField->get(i,j,k))/2;


         const real_t differential_x_0 = Nx*(-3*ConcentrationField->get(i,j,k) + 4*ConcentrationField->get(i+1,j,k) - ConcentrationField->get(i+2,j,k))/2;
         //const real_t differential_x_half = Nx*(-ConcentrationField->get(i+2,j,k) + 8*ConcentrationField->get(i+1,j,k) - 8*ConcentrationField->get(i-1,j,k) + ConcentrationField->get(i-2,j,k))/12;
         const real_t differential_x_half = Nx*(ConcentrationField->get(i+1,j,k) - ConcentrationField->get(i,j,k));
         const real_t differential_x_domain_hot = Nx*(-3*ConcentrationField->get(i,j,k) + 4*ConcentrationField->get(i+1,j,k) - ConcentrationField->get(i+2,j,k))/2;
         const real_t differential_x_domain_cold = Nx*(3*ConcentrationField->get(i,j,k) - 4*ConcentrationField->get(i-1,j,k) + ConcentrationField->get(i-2,j,k))/2;

        // const real_t differential_x_domain_hot = Nx*((real_c(-25/12))*ConcentrationField->get(i,j,k) + 4*ConcentrationField->get(i+1,j,k) - 3*ConcentrationField->get(i+2,j,k) + (real_c(4/3))*ConcentrationField->get(i+3,j,k) - (real_c(1/4))*ConcentrationField->get(i+4,j,k));
        // const real_t differential_x_domain_cold = Nx*((real_c(+25/12))*ConcentrationField->get(i,j,k) - 4*ConcentrationField->get(i+1,j,k) + 3*ConcentrationField->get(i+2,j,k) - (real_c(4/3))*ConcentrationField->get(i+3,j,k) + (real_c(1/4))*ConcentrationField->get(i+4,j,k));

         if(pos[0] == 0.5){
            nusseltx_zero += -differential_x_domain_hot;
            //WALBERLA_LOG_INFO_ON_ROOT("conc at i=0 is " << ConcentrationField->get(i,j,k));
         }

         if(pos[0] == (real_c(domainSize[0]/2) + 0.5)){

            //WALBERLA_LOG_INFO_ON_ROOT("temp at 32 is " << (ConcentrationField->get(33,j,k)) << " " << "temp at 30 is " << (ConcentrationField->get(34,j,k))  );
            nusseltx_half += ratio*VelocityField->get(*cellIt,0)*(ConcentrationField->get(*cellIt))  - differential_x_half ;

         }

         if(pos[0] == 0.5) {
            nusseltx_domain += ratio*VelocityField->get(*cellIt,0)*(ConcentrationField->get(*cellIt)) - differential_x_domain_hot;
         }
         if(pos[0] == (real_c(domainSize[0]) + 0.5))
         {
            nusseltx_domain += ratio*VelocityField->get(*cellIt,0)*(ConcentrationField->get(*cellIt)) - differential_x_domain_cold;
         }
         else{
            nusseltx_domain += ratio*VelocityField->get(*cellIt,0)*(ConcentrationField->get(*cellIt)) - differential_x_half;
         }

      }
   }
   mpi::allReduceInplace(nusseltx_domain, mpi::SUM);
   //mpi::allReduceInplace(nusseltx_zero, mpi::SUM);
   mpi::allReduceInplace(nusseltx_half, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("Nx is " << Nx << "Ny is " << Ny);
   nusseltx_domain = nusseltx_domain/(Nx*Ny*delta_theta);
   nusseltx_zero   = nusseltx_zero/(Ny*delta_theta);
   nusseltx_half   = nusseltx_half/(Nx*delta_theta);
   nusseltNumbers[0] = nusseltx_zero;
   nusseltNumbers[1] = nusseltx_half;
   nusseltNumbers[2] = nusseltx_domain;

   return nusseltNumbers;
}

real_t Simpsons (real_t bb, real_t aa, BlockDataID& ConcentrationFieldID,const shared_ptr< StructuredBlockStorage >& blocks){
   real_t nusselt_0;
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      CellInterval xyz                  = ConcentrationField->xyzSize();
      uint_t level                      = b.getLevel();
      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);
         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         const cell_idx_t i    = cellIt->x();
         const cell_idx_t j    = cellIt->y();
         const cell_idx_t k    = cellIt->z();
         if(pos[0] == real_c(0.5)){
            WALBERLA_LOG_INFO_ON_ROOT("temperature on hot wall is " << ConcentrationField->get(*cellIt));
            if(pos[1] == real_c(0.5)){
               nusselt_0 += ConcentrationField->get(*cellIt);
            }
            if(pos[1] == real_c(128.5)){
               nusselt_0 += ConcentrationField->get(*cellIt);
            }
            if(pos[1] == real_c(64.5)){
               nusselt_0 += 4 * ConcentrationField->get(*cellIt);
            }
         }
         nusselt_0 = (nusselt_0);
         //WALBERLA_LOG_INFO_ON_ROOT("nuseelt here is " << nusselt_0);

      }

   }
   return nusselt_0;
}

} // namespace walberla