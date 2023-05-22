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
//! \file NonuniformGPUCommData.impl.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/all.h"

#include "lbm_generated/gpu/NonuniformGPUCommData.h"

#include "stencil/Directions.h"

#define IDX_FLAG(d) (1 << d)

#if !defined(USE_CELL_INTERVALS)
#define INTERIOR_FLAG_BIT 29
#define INTERIOR_FLAG (1 << INTERIOR_FLAG_BIT)

#define PASS_THROUGH_FLAG_BIT 30
#define PASS_THROUGH_FLAG (1 << PASS_THROUGH_FLAG_BIT)

#define CORNER_SKIPPING_ORIGIN_FLAG_BIT 31
#define CORNER_SKIPPING_ORIGIN_FLAG (1 << CORNER_SKIPPING_ORIGIN_FLAG_BIT)
#endif

using namespace walberla::lbm_generated::util;

namespace walberla::lbm_generated {

/***********************************************************************************************************************
 *                                               Bit Mask Computation                                                  *
 **********************************************************************************************************************/

template< typename LatticeStorageSpecification_T >
void NonuniformGPUCommData< LatticeStorageSpecification_T >::registerFlags()
{
#if !defined(USE_CELL_INTERVALS)
   maskField_.registerFlag(FlagUID(true), INTERIOR_FLAG_BIT);
   maskField_.registerFlag(FlagUID(true), PASS_THROUGH_FLAG_BIT);
   maskField_.registerFlag(FlagUID(true), CORNER_SKIPPING_ORIGIN_FLAG_BIT);
#endif

   for(auto it = Stencil::beginNoCenter(); it != Stencil::end(); ++it){
      maskField_.registerFlag(FlagUID(true), Stencil::idx[*it]);
   }
}

#if defined(USE_CELL_INTERVALS)

template< typename LatticeStorageSpecification_T >
inline void NonuniformGPUCommData< LatticeStorageSpecification_T >::prepareIntervals()
{
   passThroughIntervals_.clear();
   const Block * b = dynamic_cast< const Block * >(block_);

   for(auto commDir = CommunicationStencil::beginNoCenter(); commDir != CommunicationStencil::end(); ++commDir){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*commDir);
      if(!b->neighborhoodSectionHasEquallySizedBlock(nSecIdx)){
         CellInterval ci;
         maskField_.getGhostRegion(*commDir, ci, 2);
         passThroughIntervals_.push_back(ci);
      }
   }
}

template< typename LatticeStorageSpecification_T >
inline void NonuniformGPUCommData< LatticeStorageSpecification_T >::setFlagOnInterval(const CellInterval & ci,
                                                                                   const uint_t fIdx)
{
   for(auto c : ci){
      maskField_.addFlag(c, IDX_FLAG(fIdx));
   }
}

#else

/**
 * Prepares the INTERIOR and PASS_THROUGH flags.
 * Sets the domain interior to INTERIOR. Sets any ghost layers corresponding to a coarse block
 * or no block to PASS_THROUGH.
 */
template< typename LatticeStorageSpecification_T >
void NonuniformCommData< LatticeStorageSpecification_T >::prepareFlags()
{
   const Block * b = dynamic_cast< const Block * >(block_);

   // Set interior to origin
   for (auto it = maskField_.beginXYZ(); it != maskField_.end(); ++it)
   {
      maskField_.addFlag(it.cell(), INTERIOR_FLAG);
   }

   // Set GLs to pass-through
   for(auto commDir = CommunicationStencil::beginNoCenter(); commDir != CommunicationStencil::end(); ++commDir){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*commDir);
      if(!b->neighborhoodSectionHasEquallySizedBlock(nSecIdx)){
         for(auto it = maskField_.beginGhostLayerOnlyXYZ(2, *commDir); it != maskField_.end(); ++it){
            maskField_.addFlag(it.cell(), PASS_THROUGH_FLAG);
         }
      }
   }
}

/**
 * Resets the origin flag on any ghost layers.
 */
template< typename LatticeStorageSpecification_T >
inline void NonuniformCommData< LatticeStorageSpecification_T >::resetCornerSkippingOriginFlags()
{
   const Block * b = dynamic_cast< const Block * >(block_);

   // Remove origin flag from any ghost layers
   for(auto commDir = CommunicationStencil::beginNoCenter(); commDir != CommunicationStencil::end(); ++commDir){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*commDir);
      if(!b->neighborhoodSectionHasEquallySizedBlock(nSecIdx)){
         for(auto it = maskField_.beginGhostLayerOnlyXYZ(2, *commDir); it != maskField_.end(); ++it){
            maskField_.removeFlag(it.cell(), CORNER_SKIPPING_ORIGIN_FLAG);
         }
      }
   }
}

#endif


/**
 * Determines whether the current block has the smallest BlockID among all fine blocks of a
 * given intersection volume.
 * @tparam LatticeStorageSpecification_T
 * @param cornerDir
 * @return
 */
template< typename LatticeStorageSpecification_T >
inline bool NonuniformGPUCommData< LatticeStorageSpecification_T >::haveSmallestIdInIntersection(Vector3<cell_idx_t> cornerDir)
{
   const IBlockID& myId = block_->getId();
   const Block* b = dynamic_cast< const Block* >(block_);
   return forEachSubdirectionCancel(cornerDir, [&](Vector3< cell_idx_t > dirVec) {
     const uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(dirVec[0], dirVec[1], dirVec[2]);
     if (b->neighborhoodSectionHasEquallySizedBlock(nSecIdx))
     {
        if (b->getNeighbor(nSecIdx, 0).getId() < myId) return false;
     }
     return true;
   });
}


/**
 * Sets up the feasible space for the given communication direction.
 * Additionally to the field interior, marks every ghost layer slice corresponding to an adjacent coarse block,
 * and the corresponding corner as feasible, if that corner also belongs to a coarse block and the current block
 * has the smallest BlockID participating in the intersection.
 * @param commDir A communication direction pointing toward an adjacent coarse block
 */
template< typename LatticeStorageSpecification_T >
inline void NonuniformGPUCommData< LatticeStorageSpecification_T >::setupCornerSkippingOrigins(stencil::Direction commDir)
{
#if defined(USE_CELL_INTERVALS)
   cornerSkippingOriginIntervals_.clear();
#else
   resetCornerSkippingOriginFlags();
#endif

   const Block* b = dynamic_cast< const Block* >(block_);
   Vector3<cell_idx_t> commDirVec(stencil::cx[commDir], stencil::cy[commDir], stencil::cz[commDir]);

   // Iterate all orthogonal comm directions
   forEachOrthogonalDirection< CommunicationStencil >(commDirVec, [&](Vector3< cell_idx_t > toSourceVec) {
      const uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(toSourceVec[0], toSourceVec[1], toSourceVec[2]);
      // Find if there is a coarse block or no block at all in this neighborhood
      // There are three possibilities: Coarse block, Same-level block or no block
      // Finer block is not possible because of 2:1 balance
      if (!b->neighborhoodSectionHasEquallySizedBlock(nSecIdx))
      {
         // From this adjacent coarse block (or not-block, for boundary handling), corner skipping must be handled.
         // Also, if there is no block, boundary handling in that region must be done on only
         // one of the participating fine blocks.
         Vector3< cell_idx_t > cornerDirVec = toSourceVec + commDirVec;

         // If the current block has the smallest participating ID...
         if (haveSmallestIdInIntersection(cornerDirVec))
         {
            const stencil::Direction toSourceDir = stencil::vectorToDirection(toSourceVec);

            // ... Mark source GL region as corner skipping origin.
#if defined(USE_CELL_INTERVALS)
            CellInterval ci;
            maskField_.getGhostRegion(toSourceDir, ci, 2);
            cornerSkippingOriginIntervals_.push_back(ci);
#else
            for (auto it = maskField_.beginGhostLayerOnlyXYZ(toSourceDir); it != maskField_.end(); ++it)
            {
               maskField_.addFlag(it.cell(), CORNER_SKIPPING_ORIGIN_FLAG);
            }
#endif
         }
      }
   });
}


template< typename LatticeStorageSpecification_T >
inline void NonuniformGPUCommData< LatticeStorageSpecification_T >::setupBitMaskSlice(stencil::Direction commDir, stencil::Direction streamDir)
{
   uint_t fIdx = Stencil::idx[streamDir];
   Cell streamVec(stencil::cx[streamDir], stencil::cy[streamDir], stencil::cz[streamDir]);

#if defined(USE_CELL_INTERVALS)
   CellInterval commSliceInterval;
   maskField_.getGhostRegion(commDir, commSliceInterval, 2);

   // Shift back once
   commSliceInterval.shift(-streamVec);

   // Intersect with interior and set flag on intersection volume
   CellInterval interiorIntersection(interiorInterval);
   interiorIntersection.intersect(commSliceInterval);
   if(!interiorIntersection.empty()){
      interiorIntersection.shift(streamVec);
      setFlagOnInterval(interiorIntersection, fIdx);
   }

   // Intersect with pass-through regions...
   for(auto passThroughIntersection : std::as_const(passThroughIntervals_)){
      passThroughIntersection.intersect(commSliceInterval);
      if(passThroughIntersection.empty()) continue;

      // ... shift back once more ...
      passThroughIntersection.shift(-streamVec);

      // ... intersect with interior ...
      interiorIntersection = interiorInterval;
      interiorIntersection.intersect(passThroughIntersection);
      if(!interiorIntersection.empty()){
         interiorIntersection.shift(2*streamVec.x(), 2* streamVec.y(), 2*streamVec.z());
         setFlagOnInterval(interiorIntersection, fIdx);
      }

      // ... and with corner-skipping origin regions
      for(auto originIntersection : std::as_const(cornerSkippingOriginIntervals_)){
         originIntersection.intersect(passThroughIntersection);
         if(!originIntersection.empty()){
            originIntersection.shift(2*streamVec.x(), 2* streamVec.y(), 2*streamVec.z());
            setFlagOnInterval(originIntersection, fIdx);
         }
      }
   }
#else
   for(auto it = maskField_.beginGhostLayerOnlyXYZ(2, commDir); it != maskField_.end(); ++it){
      Cell currentCell = it.cell();

      // Shift back once
      Cell shiftedCell = currentCell - streamVec;

      if (maskField_.isFlagSet(shiftedCell, INTERIOR_FLAG)){
         maskField_.addFlag(currentCell, IDX_FLAG(fIdx));
      }
      else if (maskField_.isFlagSet(shiftedCell, PASS_THROUGH_FLAG)){
         // Shift back twice
         shiftedCell -= streamVec;
         if (maskField_.isPartOfMaskSet(shiftedCell, INTERIOR_FLAG | CORNER_SKIPPING_ORIGIN_FLAG)){
            maskField_.addFlag(currentCell, IDX_FLAG(fIdx));
         }

      }
      // else continue;
   }
#endif
}

/**
 * Computes the partial coalescence bit mask on the mask field.
 * Assumes that all flags are already registered at the field, and that the field
 * has been initialized to zero.
 */
template< typename LatticeStorageSpecification_T >
void NonuniformGPUCommData< LatticeStorageSpecification_T >::computeBitMask()
{
#if defined(USE_CELL_INTERVALS)
   prepareIntervals();
#else
   prepareFlags();
#endif

   const Block* b = dynamic_cast< const Block* >(block_);
   for(auto commIt = CommunicationStencil::beginNoCenter(); commIt != CommunicationStencil::end(); ++commIt){
      stencil::Direction commDir = *commIt;
      const uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(commDir);
      if(b->neighborhoodSectionHasLargerBlock(nSecIdx)){
         setupCornerSkippingOrigins(commDir);

         for(uint_t streamDirIdx = 0; streamDirIdx < Stencil::d_per_d_length[commDir]; streamDirIdx++){
            stencil::Direction streamDir = Stencil::d_per_d[commDir][streamDirIdx];
            setupBitMaskSlice(commDir, streamDir);
         }
      }
   }
}

template< typename LatticeStorageSpecification_T >
void NonuniformGPUCommData< LatticeStorageSpecification_T >::syncDataGPU()
{
   gpu::fieldCpy(maskFieldGPU_, maskField_);
}
} // walberla::lbm_generated
