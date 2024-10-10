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
//! \file Evaluation.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "Evaluation.h"


namespace walberla
{

void Evaluation::operator()()
{
   if (checkFrequency_ == uint_c(0)) return;
   ++coarseExecutionCounter_;
   if (rampUpTime_ > coarseExecutionCounter_) return;

   if ((coarseExecutionCounter_ - uint_c(1)) % checkFrequency_ != 0) return;

   WALBERLA_CHECK_NOT_NULLPTR(blocks_)
   WALBERLA_ROOT_SECTION()
   {
      for(auto it = dragResults.begin(); it != dragResults.end(); ++it){
         coefficients_[0].push_back(it->cDRealArea);
         coefficients_[1].push_back(it->cLRealArea);
         coefficients_[2].push_back(it->cDDiscreteArea);
         coefficients_[3].push_back(it->cLDiscreteArea);
      }

      if (coefficients_[0].size() > setup_.nbrOfEvaluationPointsForCoefficientExtremas){
         for (uint_t i = uint_c(0); i < uint_c(4); ++i)
            coefficients_[i].pop_front();
      }

      for (uint_t i = uint_c(0); i < uint_c(4); ++i){
         coefficientExtremas_[i] = std::make_pair(*(coefficients_[i].begin()), *(coefficients_[i].begin()));
         for (auto v = coefficients_[i].begin(); v != coefficients_[i].end(); ++v){
            coefficientExtremas_[i].first  = std::min(coefficientExtremas_[i].first, *v);
            coefficientExtremas_[i].second = std::max(coefficientExtremas_[i].second, *v);
         }
      }

      std::ostringstream oss;

      if (logToStream_ and !dragResults.empty()){
         WALBERLA_LOG_RESULT_ON_ROOT(
            "force acting on sphere (in dimensionless lattice units of the coarsest grid - evaluated in time step "
            << coarseExecutionCounter_ - uint_c(1) << "):\n   " << force_ << oss.str()
            << "\ndrag and lift coefficients (including extrema of last " << (coefficients_[0].size() * checkFrequency_)
            << " time steps):"
               "\n   \"real\" area:"
               "\n      c_D: "
            << dragResults.back().cDRealArea << " (min = " << coefficientExtremas_[0].first << ", max = " << coefficientExtremas_[0].second
            << ")"
            << "\n      c_L: " << dragResults.back().cLRealArea << " (min = " << coefficientExtremas_[1].first
            << ", max = " << coefficientExtremas_[1].second << ")"
            << "\n   discrete area:"
               "\n      c_D: "
            << dragResults.back().cDDiscreteArea << " (min = " << coefficientExtremas_[2].first
            << ", max = " << coefficientExtremas_[2].second << ")"
            << "\n      c_L: " << dragResults.back().cLDiscreteArea << " (min = " << coefficientExtremas_[3].first
            << ", max = " << coefficientExtremas_[3].second << ")")
      }

      if (logToFile_ and !dragResults.empty()){
         std::ofstream outfile( dragFilename_.c_str(), std::ios_base::app );

         for(auto it = dragResults.begin(); it != dragResults.end(); ++it){
            outfile << it->timestep << ",";
            outfile << it->Fx << "," << it->Fy << "," << it->Fz << ",";
            outfile << it->cDRealArea << ",";
            outfile << it->cLRealArea << ",";
            outfile << it->cDDiscreteArea << ",";
            outfile << it->cLDiscreteArea;
            outfile << "\n";
         }
         outfile.close();
      }
      dragResults.clear();
   }
}

void Evaluation::resetForce()
{
   if (!initialized_) refresh();
}

void Evaluation::forceCalculation(const uint_t level)
{
   if (rampUpTime_ > coarseExecutionCounter_) return;

   if(level == maxLevel_){
      for (auto b : finestBlocks_){
         force_ += Vector3<double>(boundaryCollection_.ObstacleObject->getForce(b));
      }

      mpi::reduceInplace(force_, mpi::SUM);
      WALBERLA_ROOT_SECTION(){
         const double meanU2 = double_c(meanVelocity) * double_c(meanVelocity);

         const double cDRealArea = (double_c(8.0) * force_[0]) / (meanU2 * double_c(surfaceAreaSphere));
         const double cLRealArea = (double_c(8.0) * force_[1]) / (meanU2 * double_c(surfaceAreaSphere));
         const double cDDiscreteArea = (double_c(8.0) * force_[0]) / (meanU2 * double_c(AD_));
         const double cLDiscreteArea = (double_c(8.0) * force_[1]) / (meanU2 * double_c(AL_));

         DragCoefficient DC(fineExecutionCounter_, force_, cDRealArea, cLRealArea, cDDiscreteArea, cLDiscreteArea);
         dragResults.push_back(DC);

         fineExecutionCounter_++;
      }
   }

   force_[0] = double_c(0.0);
   force_[1] = double_c(0.0);
   force_[2] = double_c(0.0);
}

void Evaluation::refresh()
{
   WALBERLA_CHECK_NOT_NULLPTR(blocks_)
   const uint_t finestLevel = blocks_->getDepth();

   double AD(double_c(0));
   double AL(double_c(0));
   for (auto block = blocks_->begin(); block != blocks_->end(); ++block){
      const uint_t blockLevel = blocks_->getLevel(*block);
      const FlagField_T* const flagField = block->template getData< FlagField_T >(ids_.flagField);

      const auto fluid    = flagField->getFlag(setup_.fluidUID);
      const auto obstacle = flagField->getFlag(setup_.obstacleUID);
      const double area           = double_c(4.0);

      auto xyzSize = flagField->xyzSize();
      for (cell_idx_t z = xyzSize.zMin(); z <= xyzSize.zMax(); ++z){
         for (cell_idx_t y = xyzSize.yMin(); y <= xyzSize.yMax(); ++y){
            for (cell_idx_t x = xyzSize.xMin(); x <= xyzSize.xMax(); ++x){
               if (flagField->isFlagSet(x, y, z, fluid)){
                  for (auto it = Stencil_T::beginNoCenter(); it != Stencil_T::end(); ++it){
                     const cell_idx_t nx = x + cell_idx_c(it.cx());
                     const cell_idx_t ny = y + cell_idx_c(it.cy());
                     const cell_idx_t nz = z + cell_idx_c(it.cz());

                     if (flagField->isFlagSet(nx, ny, nz, obstacle)){
                        WALBERLA_CHECK(blockLevel == finestLevel, "The sphere must be completely located on the finest level")
                        if (it.cx() == 1 && it.cy() == 0 && it.cz() == 0) { AD += area; }
                        else if (it.cx() == 0 && it.cz() == 0){
                           if (it.cy() == 1){
                              AL += area;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   mpi::reduceInplace(AD, mpi::SUM);
   mpi::reduceInplace(AL, mpi::SUM);

   WALBERLA_ROOT_SECTION(){
      AD_ = AD;
      AL_ = AL;
   }
   initialized_ = true;
}
}