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
//! \file Sweeper.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include <memory>

#include "core/cell/all.h"

namespace walberla::experimental::sweep
{

class SerialSweeper
{
 private:
   std::shared_ptr< StructuredBlockForest > blocks_;

 public:
   SerialSweeper(const std::shared_ptr< StructuredBlockForest >& blocks) : blocks_{ blocks } {}

   void sweep(std::function< void(IBlock&) > func)
   {
      for (auto& block : *blocks_)
         func(block);
   }

   void sweep(std::function< void(IBlock*) > func)
   {
      for (auto& block : *blocks_)
         func(&block);
   }

   void forAllBlocks(std::function< void(IBlock&) > func) { sweep(func); }

   void forAllBlocks(std::function< void(IBlock*) > func) { sweep(func); }

   void forAllCells(std::function< void(Cell) > func)
   {
      for (cell_idx_t z = 0; z < cell_idx_c(blocks_->getNumberOfZCellsPerBlock()); ++z)
         for (cell_idx_t y = 0; y < cell_idx_c(blocks_->getNumberOfYCellsPerBlock()); ++y)
            for (cell_idx_t x = 0; x < cell_idx_c(blocks_->getNumberOfXCellsPerBlock()); ++x)
               func({ x, y, z });
   }

   void forAllCells(const CellInterval& ci, std::function< void(Cell) > func)
   {
      for (cell_idx_t z = ci.zMin(); z <= ci.zMax(); ++z)
         for (cell_idx_t y = ci.yMin(); y <= ci.yMax(); ++y)
            for (cell_idx_t x = ci.xMin(); x <= ci.xMax(); ++x)
               func({ x, y, z });
   }
};

} // namespace walberla::experimental::sweep
