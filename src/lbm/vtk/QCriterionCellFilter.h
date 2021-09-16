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
//! \file QCriterionCellFilter.h
//! \ingroup field
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/FlagField.h"
#include "core/cell/CellSet.h"
#include "core/debug/Debug.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/MacroscopicValueCalculation.h"


namespace walberla {
namespace field {

template< typename VelocityField_T, typename Filter_T >
class QCriterionCellFilter {
public:
   QCriterionCellFilter(const ConstBlockDataID velocityFieldId, Filter_T & filter, real_t lowerLimit,
         real_t upperLimit = std::numeric_limits<real_t>::infinity())
         : velocityFieldId_(velocityFieldId), filter_(filter), lowerLimit_(lowerLimit), upperLimit_(upperLimit) {}

   void operator()( CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage,
         const uint_t ghostLayers = uint_t(0) ) {
      const VelocityField_T* velocityField = block.getData< VelocityField_T >( velocityFieldId_ );
      WALBERLA_CHECK_NOT_NULLPTR(velocityField);

      const real_t dx = storage.dx(storage.getLevel(block));
      const real_t dy = storage.dy(storage.getLevel(block));
      const real_t dz = storage.dz(storage.getLevel(block));

      filter_(block);

      const cell_idx_t gl = cell_idx_c(ghostLayers);
      const cell_idx_t begin = cell_idx_c(-1) * gl;

      for (cell_idx_t z = begin; z < cell_idx_c(storage.getNumberOfZCells(block))+gl; ++z) {
         for (cell_idx_t y = begin; y < cell_idx_c(storage.getNumberOfYCells(block))+gl; ++y) {
            for (cell_idx_t x = begin; x < cell_idx_c(storage.getNumberOfXCells(block))+gl; ++x) {
               real_t q = lbm::getQCriterion(*velocityField, filter_, x, y, z, dx, dy, dz);

               if (q >= lowerLimit_ && q <= upperLimit_) {
                  filteredCells.insert(x,y,z);
               }
            }
         }
      }
   }

private:
   const ConstBlockDataID velocityFieldId_;
   Filter_T filter_;
   real_t lowerLimit_;
   real_t upperLimit_;
};

} // namespace field
} // namespace walberla
