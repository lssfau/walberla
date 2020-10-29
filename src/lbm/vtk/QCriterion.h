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
//! \file QCriterion.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/field/MacroscopicValueCalculation.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"


namespace walberla {
namespace lbm {


template< typename VelocityField_T, typename Filter_T, typename OutputType = float >
class QCriterionVTKWriter : public vtk::BlockCellDataWriter< OutputType, 1 >
{
public:
   QCriterionVTKWriter(const shared_ptr<StructuredBlockStorage> blockStorage, Filter_T & filter,
                   const ConstBlockDataID & velocityFieldId, const std::string & id ) :
         vtk::BlockCellDataWriter< OutputType, 1 >(id), blockStorage_(blockStorage), filter_(filter), velocityFieldId_(velocityFieldId), velocityField_(NULL) {}

protected:

   void configure() {
      WALBERLA_ASSERT_NOT_NULLPTR( this->block_ );
      velocityField_ = this->block_->template getData< VelocityField_T >(velocityFieldId_ );
   }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ )
   {
      WALBERLA_ASSERT_NOT_NULLPTR(velocityField_ );

      const real_t dx = blockStorage_->dx(blockStorage_->getLevel(*this->block_));
      const real_t dy = blockStorage_->dy(blockStorage_->getLevel(*this->block_));
      const real_t dz = blockStorage_->dz(blockStorage_->getLevel(*this->block_));

      filter_(*this->block_);

      return numeric_cast<OutputType>(getQCriterion(*velocityField_, filter_, x, y, z, dx, dy, dz));
   }

   const shared_ptr<StructuredBlockStorage> blockStorage_;
   Filter_T filter_;
   const ConstBlockDataID velocityFieldId_;
   const VelocityField_T* velocityField_;
};



} // namespace lbm
} // namespace walberla
