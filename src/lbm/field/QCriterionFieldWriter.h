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
//! \file VelocityFieldWriter.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/iterators/IteratorMacros.h"
#include "lbm/field/QCriterion.h"

namespace walberla {
namespace lbm {

template< typename VelocityField_T, typename QCriterionField_T, typename Filter_T >
class QCriterionFieldWriter {
public:
   QCriterionFieldWriter(const shared_ptr<StructuredBlockStorage> & blockStorage, const ConstBlockDataID & velocityFieldId,
           const BlockDataID & qCriterionFieldId, Filter_T & filter) :
      blockStorage_(blockStorage), velocityFieldId_(velocityFieldId), qCriterionFieldId_(qCriterionFieldId),
      filter_(filter)
   {}

   void operator()(IBlock * const block) {
      const VelocityField_T* velocityField = block->template getData<VelocityField_T>(velocityFieldId_);
      QCriterionField_T* qCriterionField = block->template getData<QCriterionField_T>(qCriterionFieldId_);

      WALBERLA_ASSERT_EQUAL(velocityField->xyzSize(), qCriterionField->xyzSize());

      const real_t dx = blockStorage_->dx(blockStorage_->getLevel(*block));
      const real_t dy = blockStorage_->dy(blockStorage_->getLevel(*block));
      const real_t dz = blockStorage_->dz(blockStorage_->getLevel(*block));

      filter_(*block);

      WALBERLA_FOR_ALL_CELLS_XYZ(velocityField,
         store(velocityField, qCriterionField, x, y, z, dx, dy, dz);
      )
   }

private:

   void store(const VelocityField_T* velocityField, QCriterionField_T* qCriterionField,
           const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
           const real_t dx, const real_t dy, const real_t dz) {
      qCriterionField->get(x,y,z,0) = QCriterion::get(*velocityField, filter_, x, y, z, dx, dy, dz);
   }

   const shared_ptr<StructuredBlockStorage> blockStorage_;
   ConstBlockDataID velocityFieldId_;
   BlockDataID qCriterionFieldId_;
   Filter_T & filter_;
};



} // namespace lbm
} // namespace walberla
