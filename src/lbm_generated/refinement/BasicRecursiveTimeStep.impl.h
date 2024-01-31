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
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BasicRecursiveTimeStep.h"

namespace walberla {
namespace lbm_generated {

template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T >::timestep(uint_t level)
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
      boundaryCollection_(b);
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
      boundaryCollection_(b);
      if(level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
   }

   // 2.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_){
      commScheme_->communicateFineToCoarse(level + 1);
   }
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T >::addRefinementToTimeLoop(SweepTimeloop & timeloop, uint_t level)
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

   // 2.6 Fine to Coarse Communication, receiving end
   if(level < maxLevel_)
      timeloop.addFuncBeforeTimeStep(commScheme_->communicateFineToCoarseFunctor(level + 1), "Refinement Cycle: communicate fine to coarse on level " + std::to_string(level + 1));

}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
std::function<void()> BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T >::executeStreamCollideOnLevel(uint_t level, bool withGhostLayerPropagation)
{
   return [level, withGhostLayerPropagation, this]()
   {
      if (withGhostLayerPropagation)
      {
         for(auto b: blocks_[level]){
            ghostLayerPropagation(b);
            sweepCollection_.streamCollide(b);
         }
      }
      else
      {
         for(auto b: blocks_[level]){
            sweepCollection_.streamCollide(b);
         }
      }
   };
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
std::function<void()>  BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T >::executeBoundaryHandlingOnLevel(uint_t level)
{
   return [level, this]() {
      for (auto b : blocks_[level])
      {
         boundaryCollection_(b);
         if (level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
      }
   };
}


template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T >
void BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T >::ghostLayerPropagation(Block * block)
{
   auto pdfField = block->getData<PdfField_T>(pdfFieldId_);

   for(auto it = CommunicationStencil::beginNoCenter(); it != CommunicationStencil::end(); ++it){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*it);
      // Propagate on ghost layers shadowing coarse or no blocks
      if(block->neighborhoodSectionHasLargerBlock(nSecIdx)){
         CellInterval ci;
         pdfField->getGhostRegion(*it, ci, 1);
         sweepCollection_.streamOnlyNoAdvancementCellInterval(block, ci);
      }
   }
}

// Refinement Timestep from post collision state:
//template< typename PdfField_T, typename LbSweep_T >
//void BasicRecursiveTimeStep< PdfField_T, LbSweep_T >::timestep(uint_t level)
//{
//   std::vector<Block *> blocks;
//   sbfs_->getBlocks(blocks, level);
//
//   uint_t maxLevel = sbfs_->getDepth();
//
//   // 1.1 Equal-Level Communication
//   commScheme_->communicateEqualLevel(level);
//
//   // 1.2 Coarse to Fine Communication
//   if(level < maxLevel){
//      commScheme_->communicateCoarseToFine(level + 1);
//   }
//
//   // 1.3 Boundary Handling and
//   // 1.4 Prepare Coalescence (which happens during the recursive descent)
//   for(auto b : blocks){
//      boundaryFunctor_(b);
//      if(level != maxLevel) pdfFieldPackInfo_->prepareCoalescence(b);
//   }
//
//   // 1.5 Recursive Descent
//   if(level < maxLevel){
//      timestep(level + 1);
//   }
//
//   // 1.6 First Collision and ghost-layer propagation
//   for(auto b: blocks){
//      if(level != 0) ghostLayerPropagation(b);  // GL-Propagation first without swapping arrays...
//      sweepCollection_.streamCollide(b);                // then Stream-Collide on interior, and swap arrays
//   }
//
//   // Stop here if on coarsest level.
//   // Otherwise, continue to second subcycle.
//   if(level == 0) return;
//
//   // 2.1 Equal-Level Communication
//   commScheme_->communicateEqualLevel(level);
//
//   // 2.2 Coarse to Fine Communication
//   if(level < maxLevel){
//      commScheme_->communicateCoarseToFine(level + 1);
//   }
//
//   // 2.3 Boundary Handling and
//   // 2.4 Prepare Coalescence (which happens during the recursive descent)
//   for(auto b : blocks){
//      boundaryFunctor_(b);
//      if(level != maxLevel) pdfFieldPackInfo_->prepareCoalescence(b);
//   }
//
//   // 2.5 Recursive Descent
//   if(level < maxLevel){
//      timestep(level + 1);
//   }
//
//   // 2.6 Fine to Coarse Communication
//   commScheme_->communicateFineToCoarse(level);
//
//   // 2.7 Second Collision
//   for(auto b: blocks){
//      sweepCollection_.streamCollide(b);
//   }
//}

} // namespace lbm_generated
} // namespace walberla
