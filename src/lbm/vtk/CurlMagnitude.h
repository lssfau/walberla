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
//! \file Velocity.h
//! \ingroup lbm
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"


namespace walberla {
namespace lbm {

/**
 * lengthScaleWeight is only used to be able to derive quantities for the adaptive mesh refinement criterium of curl.
 *
 * @tparam VelocityField_T
 * @tparam Filter_T
 * @tparam OutputType
 */
template< typename VelocityField_T, typename Filter_T, typename OutputType = float >
class CurlMagnitudeVTKWriter : public vtk::BlockCellDataWriter< OutputType, 1 >
{
public:
   CurlMagnitudeVTKWriter(const shared_ptr<StructuredBlockStorage> blockStorage, Filter_T & filter,
         const ConstBlockDataID & velocityFieldId, const std::string & id, const real_t lengthScaleWeight = real_t(-1)) :
         vtk::BlockCellDataWriter< OutputType, 1 >(id), blockStorage_(blockStorage), filter_(filter),
         velocityFieldId_(velocityFieldId), velocityField_(NULL), lengthScaleWeight_(lengthScaleWeight) {}

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

      const auto one = cell_idx_t(1);

      const auto halfInvDx = real_t(0.5) * real_t(1) / dx;
      const auto halfInvDy = real_t(0.5) * real_t(1) / dy;
      const auto halfInvDz = real_t(0.5) * real_t(1) / dz;

      const real_t lengthScale = std::cbrt(dx*dy*dz);
      const real_t weightedLengthScale = !realIsIdentical(lengthScaleWeight_, real_t(-1)) ?
            std::pow(lengthScale, (lengthScaleWeight_+1)/lengthScaleWeight_) : real_t(1);

      if(filter_(x,y,z) && filter_(x+one,y,z) && filter_(x-one,y,z) && filter_(x,y+one,z)
         && filter_(x,y-one,z) && filter_(x,y,z+one) && filter_(x,y,z-one)) {
         const Vector3< real_t > xa = velocityField_->get(x+one,y,z);
         const Vector3< real_t > xb = velocityField_->get(x-one,y,z);
         const Vector3< real_t > ya = velocityField_->get(x,y+one,z);
         const Vector3< real_t > yb = velocityField_->get(x,y-one,z);
         const Vector3< real_t > za = velocityField_->get(x,y,z+one);
         const Vector3< real_t > zb = velocityField_->get(x,y,z-one);

         const real_t duxdy = (ya[0] - yb[0]) * halfInvDy;
         const real_t duxdz = (za[0] - zb[0]) * halfInvDz;

         const real_t duydx = (xa[1] - xb[1]) * halfInvDx;
         const real_t duydz = (za[1] - zb[1]) * halfInvDz;

         const real_t duzdx = (xa[2] - xb[2]) * halfInvDx;
         const real_t duzdy = (ya[2] - yb[2]) * halfInvDy;

         const Vector3< real_t > curl( duzdy - duydz, duxdz - duzdx, duydx - duxdy );
         const auto curlMag = curl.length();

         const auto curlSensor = curlMag * weightedLengthScale;

         return numeric_cast< OutputType >(curlSensor);
      }

      return numeric_cast< OutputType >(0);
   }

   const shared_ptr<StructuredBlockStorage> blockStorage_;
   Filter_T filter_;
   const ConstBlockDataID velocityFieldId_;
   const VelocityField_T* velocityField_;
   const real_t lengthScaleWeight_;
};



} // namespace lbm
} // namespace walberla
