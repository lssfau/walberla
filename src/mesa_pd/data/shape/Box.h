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
//! \file Box.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class Box : public BaseShape
{
public:
   explicit Box(const Vec3& edgeLength)
      : BaseShape(Box::SHAPE_TYPE)
      , edgeLength_(edgeLength)
   {}

   const Vec3&   getEdgeLength() const   { return edgeLength_; }

   real_t getVolume() const override { return edgeLength_[0] * edgeLength_[1] * edgeLength_[2]; };
   void   updateMassAndInertia(const real_t density) override;

   static const int SHAPE_TYPE = 3; ///< Unique shape type identifier for boxes.\ingroup mesa_pd_shape

private:
   Vec3   edgeLength_;   ///< edge length of the box
};

inline
void Box::updateMassAndInertia(const real_t density)
{
   const real_t m = getVolume() * density;
   const Mat3   I = Mat3::makeDiagonalMatrix(
         edgeLength_[1]*edgeLength_[1] + edgeLength_[2]*edgeLength_[2] ,
         edgeLength_[0]*edgeLength_[0] + edgeLength_[2]*edgeLength_[2] ,
         edgeLength_[0]*edgeLength_[0] + edgeLength_[1]*edgeLength_[1] ) * (m / static_cast<real_t>( 12 ));

   getInvMass()      = real_t(1.0) / m;
   getInvInertiaBF() = I.getInverse();
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
