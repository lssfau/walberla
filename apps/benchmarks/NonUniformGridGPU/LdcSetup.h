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
//! \file LdcSetup.h
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================
#pragma once

#include "blockforest/SetupBlock.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/all.h"

#include "field/FlagField.h"

#include "field/FlagUID.h"

using namespace walberla;

using RefinementSelectionFunctor = SetupBlockForest::RefinementSelectionFunction;

using FlagField_T          = FlagField< uint8_t >;

class LDCRefinement
{
 private:
   const uint_t refinementDepth_;

 public:
   explicit LDCRefinement(const uint_t depth) : refinementDepth_(depth){};

   void operator()(SetupBlockForest& forest) const
   {
      const AABB & domain = forest.getDomain();

      const AABB leftCorner( 0, domain.yMax() -1, 0, 1, domain.yMax() , domain.zMax() );
      const AABB rightCorner( domain.xMax() - 1, domain.yMax() -1, 0, domain.xMax(), domain.yMax() , domain.zMax() );

      for(auto & block : forest)
      {
         auto & aabb = block.getAABB();
         if( leftCorner.intersects( aabb ) || rightCorner.intersects( aabb ) )
         {
            if( block.getLevel() < refinementDepth_)
               block.setMarker( true );
         }
      }
   }
};

class LDC
{
 private:
   const std::string refinementProfile_;
   const uint_t refinementDepth_;

   const FlagUID noSlipFlagUID_;
   const FlagUID ubbFlagUID_;

 public:
   explicit LDC(const uint_t depth) : refinementDepth_(depth), noSlipFlagUID_("NoSlip"), ubbFlagUID_("UBB"){};

   RefinementSelectionFunctor refinementSelector() const
   {
      return LDCRefinement(refinementDepth_);
   }

   void setupBoundaryFlagField(StructuredBlockForest& sbfs, const BlockDataID flagFieldID)
   {
      for (auto bIt = sbfs.begin(); bIt != sbfs.end(); ++bIt)
      {
         auto& b           = dynamic_cast< Block& >(*bIt);
         const uint_t level       = b.getLevel();
         auto flagField     = b.getData< FlagField_T >(flagFieldID);
         const uint8_t noslipFlag = flagField->registerFlag(noSlipFlagUID_);
         const uint8_t ubbFlag    = flagField->registerFlag(ubbFlagUID_);
         for (auto cIt = flagField->beginWithGhostLayerXYZ(2); cIt != flagField->end(); ++cIt)
         {
            const Cell localCell = cIt.cell();
            Cell globalCell(localCell);
            sbfs.transformBlockLocalToGlobalCell(globalCell, b);
            if (globalCell.y() >= cell_idx_c(sbfs.getNumberOfYCells(level))) { flagField->addFlag(localCell, ubbFlag); }
            else if (globalCell.y() < 0 || globalCell.x() < 0 || globalCell.x() >= cell_idx_c(sbfs.getNumberOfXCells(level)))
            {
               flagField->addFlag(localCell, noslipFlag);
            }
         }
      }
   }
};
