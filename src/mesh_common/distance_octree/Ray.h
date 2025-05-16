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
//! \file Ray.h
//! \ingroup mesh
//! \author Lukas Werner
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#include "core/math/Vector3.h"

namespace walberla {
namespace mesh {

class Ray {
 public:

   Ray () : Ray (Vector3<real_t>(0,0,0), Vector3<real_t>(1,0,0)) {}

   Ray (Vector3<real_t> origin, Vector3<real_t> direction) {
      setDirection(direction);
      setOrigin(origin);
   }

 private:

   Vector3<real_t> origin_; //!< Origin of the ray.
   Vector3<real_t> direction_; //!< The normalized direction of the ray.
   Vector3<real_t> inv_direction_; //!< The inverted direction of the ray.
   Vector3<int8_t> sign_; /*!< The signs of the inverted direction of the ray.
                             (Required for Ray-Box intersection code.)*/
   size_t imageX_; //!< Y value of the pixel coordinate this ray intersects.
   size_t imageY_; //!< X value of the pixel coordinate this ray intersects.

 public:

   inline const Vector3<real_t>& getOrigin () const {
      return origin_;
   }

   inline const Vector3<real_t>& getDirection () const {
      return direction_;
   }

   inline real_t getDirection (size_t axis) const {
      WALBERLA_ASSERT(axis <= 2, "No valid axis index passed.");
      return direction_[axis];
   }

   inline real_t xDir () const {
      return direction_[0];
   }

   inline real_t yDir () const {
      return direction_[1];
   }

   inline real_t zDir () const {
      return direction_[2];
   }

   inline const Vector3<real_t>& getInvDirection () const {
      return inv_direction_;
   }

   inline real_t getInvDirection (size_t axis) const {
      WALBERLA_ASSERT(axis <= 2, "No valid axis index passed.");
      return inv_direction_[axis];
   }

   inline real_t xInvDir () const {
      return inv_direction_[0];
   }

   inline real_t yInvDir () const {
      return inv_direction_[1];
   }

   inline real_t zInvDir () const {
      return inv_direction_[2];
   }

   inline const Vector3<int8_t>& getInvDirectionSigns () const {
      return sign_;
   }

   inline size_t getImageX () const {
      return imageX_;
   }

   inline size_t getImageY () const {
      return imageY_;
   }

   inline void setOrigin (const Vector3<real_t>& origin) {
      origin_ = origin;
   }

   inline void setDirection (const Vector3<real_t>& direction) {
      WALBERLA_CHECK_FLOAT_EQUAL(direction.length(), real_t(1));
      direction_ = direction;
      calcInvDirection();
   }

   inline void setImageCoordinate (size_t x, size_t y) {
      imageX_ = x;
      imageY_ = y;
   }

   inline void setImageX (size_t x) {
      imageX_ = x;
   }

   inline void setImageY (size_t y) {
      imageY_ = y;
   }

   inline void calcInvDirection () {
      inv_direction_ = Vector3<real_t>(1/direction_[0], 1/direction_[1], 1/direction_[2]);
      sign_[0] = (inv_direction_[0] < 0) ? int8_t(1) : int8_t(0);
      sign_[1] = (inv_direction_[1] < 0) ? int8_t(1) : int8_t(0);
      sign_[2] = (inv_direction_[2] < 0) ? int8_t(1) : int8_t(0);
   }
};


inline bool rayAABBIntersection(const AABB& aabb, const Ray& ray, real_t& t, real_t padding, Vector3<real_t>* n) {
   // An Efficient and Robust Ray–Box Intersection Algorithm: http://people.csail.mit.edu/amy/papers/box-jgt.pdf
   const Vector3<real_t> paddingVector(padding, padding, padding);
   Vector3<real_t> bounds[2] = {
      aabb.min() - paddingVector,
      aabb.max() + paddingVector
   };

   const Vector3<int8_t>& sign = ray.getInvDirectionSigns();
   const Vector3<real_t>& invDirection = ray.getInvDirection();
   const Vector3<real_t>& origin = ray.getOrigin();

   const real_t realMax = std::numeric_limits<real_t>::max();

   size_t tminAxis = 0, tmaxAxis = 0;
   real_t txmin, txmax;
   real_t tmin = txmin = (bounds[sign[0]][0] - origin[0]) * invDirection[0];
   real_t tmax = txmax = (bounds[1-sign[0]][0] - origin[0]) * invDirection[0];
   real_t tymin = (bounds[sign[1]][1] - origin[1]) * invDirection[1];
   real_t tymax = (bounds[1-sign[1]][1] - origin[1]) * invDirection[1];
   if (tmin > tymax || tymin > tmax) {
      t = realMax;
      return false;
   }
   if (tymin > tmin) {
      tminAxis = 1;
      tmin = tymin;
   }
   if (tymax < tmax) {
      tmaxAxis = 1;
      tmax = tymax;
   }
   real_t tzmin = (bounds[sign[2]][2] - origin[2]) * invDirection[2];
   real_t tzmax = (bounds[1-sign[2]][2] - origin[2]) * invDirection[2];
   if (tmin > tzmax || tzmin > tmax) {
      t = realMax;
      return false;
   }
   if (tzmin > tmin) {
      tminAxis = 2;
      tmin = tzmin;
   }
   if (tzmax < tmax) {
      tmaxAxis = 2;
      tmax = tzmax;
   }

   if (n != nullptr) {
      (*n)[0] = (*n)[1] = (*n)[2] = real_t(0);
   }
   real_t t_;
   if (tmin > 0) {
      // ray hit box from outside
      t_ = tmin;
      if (n != nullptr) {
         (*n)[tminAxis] = real_t(1);
      }
   } else if (tmax < 0) {
      // tmin and tmax are smaller than 0 -> box is in rays negative direction
      t = realMax;
      return false;
   } else {
      // ray origin within box
      t_ = tmax;
      if (n != nullptr) {
         (*n)[tmaxAxis] = real_t(1);
      }
   }

   if (n != nullptr) {
      if (ray.getDirection() * (*n) > 0) {
         *n = -(*n);
      }
   }

   t = t_;
   return true;
}



} //namespace mesh
} //namespace walberla
