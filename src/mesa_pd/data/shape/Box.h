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
   explicit Box(const Vec3& edgeLength = Vec3(real_t(1)))
      : BaseShape(Box::SHAPE_TYPE)
      , edgeLength_(edgeLength)
   {}

   const Vec3&   getEdgeLength() const   { return edgeLength_; }

   real_t getVolume() const override { return edgeLength_[0] * edgeLength_[1] * edgeLength_[2]; };
   void   updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& bfD ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 3; ///< Unique shape type identifier for boxes.\ingroup mesa_pd_shape

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

   mass_ = m;
   invMass_      = real_t(1.0) / m;

   inertiaBF_ = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Box::support( const Vec3& bfD ) const
{
   //As it is save to say we have atleast one component of the d-vector != 0 we can use
   Vec3 relativSupport = Vec3( math::sign(bfD[0])*edgeLength_[0]*real_t(0.5),
                               math::sign(bfD[1])*edgeLength_[1]*real_t(0.5),
                               math::sign(bfD[2])*edgeLength_[2]*real_t(0.5) );

   return relativSupport;
}

inline
void Box::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << edgeLength_;
}
inline
void Box::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> edgeLength_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
