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
//! \file UpdateSecondGhostLayer.h
//! \ingroup blockforest
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

namespace walberla
{
namespace blockforest
{
/********************************************************************************************************************
 * Manual update of a field's second and further ghost layers in setups with domain size = 1 and periodicity in at
 * least one direction.
 *
 * If a field has two ghost layers and if the inner field has only a size of one (in one or more directions),
 * regular communication does not update the second (and further) ghost layers correctly. Here, the content of the
 * first ghost layers is manually copied to the second (and further) ghost layers when the dimension of the field is
 * one in this direction.
 *******************************************************************************************************************/
template< typename Field_T >
class UpdateSecondGhostLayer
{
 public:
   UpdateSecondGhostLayer(const std::weak_ptr< StructuredBlockForest >& blockForestPtr, BlockDataID fieldID)
      : blockForest_(blockForestPtr), fieldID_(fieldID)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         Field_T* const field = blockIt->template getData< Field_T >(fieldID_);

         // this function is only necessary if the flag field has at least two ghost layers
         if (field->nrOfGhostLayers() < uint_c(2)) { isNecessary_ = false; }
         else
         {
            // running this function is only necessary when the field is of size 1 in at least one direction
            isNecessary_ = (field->xSize() == uint_c(1) && blockForest->isXPeriodic()) ||
                           (field->ySize() == uint_c(1) && blockForest->isYPeriodic()) ||
                           (field->zSize() == uint_c(1) && blockForest->isZPeriodic());
         }
      }
   }

   void operator()()
   {
      if (!isNecessary_) { return; }

      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         Field_T* const field = blockIt->template getData< Field_T >(fieldID_);

         if (field->xSize() == uint_c(1) && blockForest->isXPeriodic())
         {
            // iterate ghost layer at x == -1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::W, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at x == -1 to x == -2
               fieldIt.neighbor(stencil::W) = *fieldIt;
            }

            // iterate ghost layer at x == xSize()+1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::E, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at x == xSize()+1 to x == xSize()+2
               fieldIt.neighbor(stencil::E) = *fieldIt;
            }
         }

         if (field->ySize() == uint_c(1) && blockForest->isYPeriodic())
         {
            // iterate ghost layer at y == -1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::S, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at y == -1 to y == -2
               fieldIt.neighbor(stencil::S) = *fieldIt;
            }

            // iterate ghost layer at y == ySize()+1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::N, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at y == ySize()+1 to y == ySize()+2
               fieldIt.neighbor(stencil::N) = *fieldIt;
            }
         }

         if (field->zSize() == uint_c(1) && blockForest->isZPeriodic())
         {
            // iterate ghost layer at z == -1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::B, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at z == -1 to z == -2
               fieldIt.neighbor(stencil::B) = *fieldIt;
            }

            // iterate ghost layer at y == zSize()+1
            for (auto fieldIt = field->beginGhostLayerOnly(uint_c(1), stencil::T, true); fieldIt != field->end();
                 ++fieldIt)
            {
               // copy data from ghost layer at z == zSize()+1 to z == zSize()+2
               fieldIt.neighbor(stencil::T) = *fieldIt;
            }
         }
      }
   }

 private:
   std::weak_ptr< StructuredBlockForest > blockForest_;

   BlockDataID fieldID_;

   bool isNecessary_;
}; // class UpdateSecondGhostLayer

} // namespace blockforest
} // namespace walberla
