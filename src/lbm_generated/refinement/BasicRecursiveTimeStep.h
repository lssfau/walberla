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
//! \file BasicRecursiveTimeStep.h
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/NonUniformBufferedScheme.h"

#include "lbm/field/PdfField.h"
#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla {

using blockforest::communication::NonUniformBufferedScheme;
using BlockFunction = std::function<void (const uint_t)>; // parameters: block, level

namespace lbm_generated {

/**
 *
 * @tparam LatticeStorageSpecification_T   Generated storage specification
 * @tparam SweepCollection_T LBM SweepCollection (must be able to call stream, collide, streamCollide and streamOnlyNoAdvancement)
 * @tparam BoundaryCollection_T LBM Boundary collection (Functor that runs all boundary kernels at call)
 */
template< typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
class BasicRecursiveTimeStep
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using Stencil = typename LatticeStorageSpecification_T::Stencil;
   using CommunicationStencil = typename LatticeStorageSpecification_T::CommunicationStencil;
   using CommScheme = NonUniformBufferedScheme< CommunicationStencil >;
   using PackInfo = lbm_generated::NonuniformGeneratedPdfPackInfo< PdfField_T >;

   BasicRecursiveTimeStep(std::shared_ptr< StructuredBlockForest > & sbfs,
                          const BlockDataID & pdfFieldId, SweepCollection_T & sweepCollection, BoundaryCollection_T & boundaryCollection,
                          std::shared_ptr< CommScheme > & commScheme, std::shared_ptr< PackInfo > & pdfFieldPackInfo):
      sbfs_(sbfs), pdfFieldId_(pdfFieldId), pdfFieldPackInfo_(pdfFieldPackInfo), commScheme_(commScheme),
      sweepCollection_(sweepCollection), boundaryCollection_(boundaryCollection)
      {
#ifndef NDEBUG
      for (auto& block : *sbfs)
         WALBERLA_ASSERT(block.isDataOfType<PdfField_T>(pdfFieldId_), "Template parameter PdfField_T is of different type than BlockDataID pdfFieldId that is provided as constructor argument")
#endif
      maxLevel_ = sbfs->getDepth();

      for (uint_t level = 0; level <= maxLevel_; level++)
      {
         std::vector<Block *> blocks;
         sbfs->getBlocks(blocks, level);
         blocks_.push_back(blocks);
      }
     };

   void operator() () { timestep(0); };
   void addRefinementToTimeLoop(SweepTimeloop & timeloop, uint_t level=0);

   inline void addPostBoundaryHandlingBlockFunction( const BlockFunction & function );

 private:
   void timestep(uint_t level);
   void ghostLayerPropagation(Block * block);
   std::function< void() > executeStreamCollideOnLevel(uint_t level, bool withGhostLayerPropagation=false);
   std::function< void() > executeBoundaryHandlingOnLevel(uint_t level);
   std::function< void() > executePostBoundaryBlockFunctions(uint_t level);

   std::shared_ptr< StructuredBlockForest > sbfs_;
   uint_t maxLevel_;
   std::vector<std::vector<Block *>> blocks_;

   const BlockDataID pdfFieldId_;
   std::shared_ptr< PackInfo > pdfFieldPackInfo_;
   std::shared_ptr< CommScheme > commScheme_;

   SweepCollection_T & sweepCollection_;
   BoundaryCollection_T & boundaryCollection_;

   std::vector< BlockFunction >  globalPostBoundaryHandlingBlockFunctions_;
};

} // namespace lbm_generated
} // namespace walberla

#include "lbm_generated/refinement/BasicRecursiveTimeStep.impl.h"
