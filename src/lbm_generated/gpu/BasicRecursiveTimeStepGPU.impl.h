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
//! \file BasicRecursiveTimeStep.impl.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BasicRecursiveTimeStepGPU.h"

namespace walberla {
namespace lbm_generated {

template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::timestep(uint_t level)
{
   // 1.1 Collision
   for(auto b: blocks_[level]){
      sweepCollection_.streamCollide(b);
   }

   // 1.2 Recursive Descent
   if(level < maxLevel_){
      timestep(level + 1);
   }

   // 1.3 Coarse to Fine Communication, receiving end
   if(level != 0){
      commScheme_->communicateCoarseToFine(level);
   }

   // 1.4 Equal-Level Communication
   commScheme_->communicateEqualLevel(level);

   // 1.5 Boundary Handling and Coalescence Preparation
   for(auto b : blocks_[level]){
      boundaryCollection_(b, nullptr);
      if(level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
   }

   // 1.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_){
      commScheme_->communicateFineToCoarse(level + 1);
   }

   // Stop here if on coarsest level.
   // Otherwise, continue to second subcycle.
   if(level == 0) return;

   // 2.1 Collision and Ghost-Layer Propagation
   for(auto b: blocks_[level]){
      ghostLayerPropagation(b);  // GL-Propagation first without swapping arrays...
      sweepCollection_.streamCollide(b);                // then Stream-Collide on interior, and swap arrays
   }

   // 2.2 Recursive Descent
   if(level < maxLevel_){
      timestep(level + 1);
   }

   // 2.4 Equal-Level Communication
   commScheme_->communicateEqualLevel(level);

   // 2.5 Boundary Handling and Coalescence Preparation
   for(auto b : blocks_[level]){
      boundaryCollection_(b, nullptr);
      if(level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
   }

   // 2.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_){
      commScheme_->communicateFineToCoarse(level + 1);
   }
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::addRefinementToTimeLoop(SweepTimeloop & timeloop, uint_t level)
{
   // 1.1 Collision
   timeloop.addFuncBeforeTimeStep(executeStreamCollideOnLevel(level), "Refinement Cycle: streamCollide on level " + std::to_string(level));

   // 1.2 Recursive Descent
   if(level < maxLevel_){
      addRefinementToTimeLoop(timeloop, level + 1);
   }

   // 1.3 Coarse to Fine Communication, receiving end
   if(level != 0){
      timeloop.addFuncBeforeTimeStep(commScheme_->communicateCoarseToFineFunctor(level), "Refinement Cycle: communicate coarse to fine on level " + std::to_string(level));
   }

   // 1.4 Equal-Level Communication
   timeloop.addFuncBeforeTimeStep(commScheme_->communicateEqualLevelFunctor(level), "Refinement Cycle: communicate equal level on level " + std::to_string(level));


   // 1.5 Boundary Handling and Coalescence Preparation
   timeloop.addFuncBeforeTimeStep(executeBoundaryHandlingOnLevel(level), "Refinement Cycle: boundary handling on level " + std::to_string(level));
   timeloop.addFuncBeforeTimeStep(executePostBoundaryBlockFunctions(level), "Refinement Cycle: post boundary handling block functions on level " + std::to_string(level));

   // 1.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_){
      timeloop.addFuncBeforeTimeStep(commScheme_->communicateFineToCoarseFunctor(level + 1), "Refinement Cycle: communicate fine to coarse on level " + std::to_string(level + 1));
   }

   // Stop here if on coarsest level.
   // Otherwise, continue to second subcycle.
   if(level == 0) return;

   // 2.1 Collision and Ghost-Layer Propagation
   timeloop.addFuncBeforeTimeStep(executeStreamCollideOnLevel(level, true), "Refinement Cycle: streamCollide with ghost layer propagation on level " + std::to_string(level));

   // 2.2 Recursive Descent
   if(level < maxLevel_)
      addRefinementToTimeLoop(timeloop, level + 1);


   // 2.4 Equal-Level Communication
   timeloop.addFuncBeforeTimeStep(commScheme_->communicateEqualLevelFunctor(level), "Refinement Cycle: communicate equal level on level " + std::to_string(level));

   // 2.5 Boundary Handling and Coalescence Preparation
   timeloop.addFuncBeforeTimeStep(executeBoundaryHandlingOnLevel(level), "Refinement Cycle: boundary handling on level " + std::to_string(level));
   timeloop.addFuncBeforeTimeStep(executePostBoundaryBlockFunctions(level), "Refinement Cycle: post boundary handling block functions on level " + std::to_string(level));

   // 2.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_)
      timeloop.addFuncBeforeTimeStep(commScheme_->communicateFineToCoarseFunctor(level + 1), "Refinement Cycle: communicate fine to coarse on level " + std::to_string(level + 1));

}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
std::function<void()> BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::executeStreamCollideOnLevel(uint_t level, bool withGhostLayerPropagation)
{
   if(sweepCollection_.blockWise()){
      return [level, withGhostLayerPropagation, this](){
         if (withGhostLayerPropagation){
            const uint8_t timestepPlusOne = (timestepPerLevel_[level] + 1) & 1;
            sweepCollection_.ghostLayerPropagation(level, timestepPlusOne);
            WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
            WALBERLA_GPU_CHECK(gpuPeekAtLastError())
            timestepPerLevel_[level] = (timestepPerLevel_[level] + 1) & 1;
            sweepCollection_.streamCollideOverBlocks(level, timestepPerLevel_[level]);
            for (uint_t i = 0; i < blocks_[level].size(); i++){
               auto pdfs = blocks_[level][i]->getData< PdfField_T >(pdfFieldId_);
               pdfs->advanceTimestep();
            }
            commScheme_->setTimestepForLevel(level, timestepPerLevel_[level]);
            WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
            WALBERLA_GPU_CHECK(gpuPeekAtLastError())
         }
         else{
            timestepPerLevel_[level] = (timestepPerLevel_[level] + 1) & 1;
            sweepCollection_.streamCollideOverBlocks(level, timestepPerLevel_[level]);
            for (uint_t i = 0; i < blocks_[level].size(); i++){
               auto pdfs = blocks_[level][i]->getData< PdfField_T >(pdfFieldId_);
               pdfs->advanceTimestep();
            }
            commScheme_->setTimestepForLevel(level, timestepPerLevel_[level]);
            WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
            WALBERLA_GPU_CHECK(gpuPeekAtLastError())
         }
      };
   }
   else{
      return [level, withGhostLayerPropagation, this](){
         if (withGhostLayerPropagation){
            for (uint_t i = 0; i < blocks_[level].size(); i++){
               ghostLayerPropagation(blocks_[level][i], streams_[level][i % nStreams_]);
               sweepCollection_.streamCollide(blocks_[level][i], 0, streams_[level][i % nStreams_]);
            }
         }
         else{
            for (uint_t i = 0; i < blocks_[level].size(); i++){
               sweepCollection_.streamCollide(blocks_[level][i], 0, streams_[level][i % nStreams_]);
            }
         }
      };
   }
}

template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
std::function<void()> BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::executeBoundaryHandlingOnLevel(uint_t level)
{
   return [this, level]() {
      for (uint_t i = 0; i < blocks_[level].size(); i++){
         boundaryCollection_(blocks_[level][i], streams_[level][i % nStreams_]);
         if (level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(blocks_[level][i], streams_[level][i % nStreams_]);
      }
   };
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::ghostLayerPropagation(
   Block * block, gpuStream_t gpuStream)
{
   auto pdfField = block->getData<PdfField_T>(pdfFieldId_);

   for(auto it = CommunicationStencil::beginNoCenter(); it != CommunicationStencil::end(); ++it){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*it);
      // Propagate on ghost layers shadowing coarse or no blocks
      if(block->neighborhoodSectionHasLargerBlock(nSecIdx)){
         CellInterval ci;
         pdfField->getGhostRegion(*it, ci, 1);
         sweepCollection_.streamOnlyNoAdvancementCellInterval(block, ci, gpuStream);
      }
   }
}

template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
std::function<void()> BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::executePostBoundaryBlockFunctions(uint_t level)
{
   return [this, level]() {
      for( const auto& func : globalPostBoundaryHandlingBlockFunctions_ ){
         func(level);
      }
   };
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
inline void BasicRecursiveTimeStepGPU< PdfField_T, SweepCollection_T, BoundaryCollection_T >::addPostBoundaryHandlingBlockFunction( const BlockFunction & function )
{
   globalPostBoundaryHandlingBlockFunctions_.emplace_back( function );
}


} // namespace lbm_generated
} // namespace walberla
