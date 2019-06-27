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
//! \file BaseShape.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * Base class of all shape types.
 *
 * This class contains the basic properties which are inherent for all shapes.
 * It also stores the shapeType identifier which allows to cast to the actual shape.
 */
class BaseShape
{
public:
   explicit BaseShape(const int shapeType) : shapeType_(shapeType) {}
   virtual ~BaseShape() = default;

   ///Updates mass and inertia according to the actual shape.
   virtual void updateMassAndInertia(const real_t density);

   virtual real_t getVolume() const = 0;

   const real_t& getMass() const {return mass_;}
   const real_t& getInvMass() const {return invMass_;}

   const Mat3& getInertiaBF() const {return inertiaBF_;}
   const Mat3& getInvInertiaBF() const {return invInertiaBF_;}

   const int& getShapeType() const {return shapeType_;}

   virtual Vec3 support( const Vec3& /*d*/ ) const {WALBERLA_ABORT("Not implemented!");}

   static const int INVALID_SHAPE = -1; ///< Unique *invalid* shape type identifier.\ingroup mesa_pd_shape
protected:
   real_t mass_         = real_t(0);       ///< mass
   real_t invMass_      = real_t(0);       ///< inverse mass
   Mat3   inertiaBF_    = Mat3(real_t(0)); ///< inertia matrix in the body frame
   Mat3   invInertiaBF_ = Mat3(real_t(0)); ///< inverse inertia matrix in the body frame
   int    shapeType_    = INVALID_SHAPE;   ///< \ingroup mesa_pd_shape
};

inline
void BaseShape::updateMassAndInertia(const real_t /*density*/)
{
   WALBERLA_ABORT("updateMassAndInertia not implemented!");
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
