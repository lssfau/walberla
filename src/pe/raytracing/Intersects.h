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
//! \file Intersects.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/Types.h>
#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/CylindricalBoundary.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"
#include "pe/utility/BodyCast.h"
#include <boost/math/special_functions/sign.hpp>
#include <boost/tuple/tuple.hpp>

#include <pe/raytracing/Ray.h>

#define EPSILON real_t(1e-4)

using namespace boost::math;

namespace walberla {
namespace pe {
namespace raytracing {
inline bool intersects(const AABB& aabb, const Ray& ray, real_t& t, real_t padding = real_t(0.0));

inline bool intersects(const SphereID sphere, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const PlaneID plane, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const BoxID box, const Ray& ray, real_t& t, Vec3& n);

struct IntersectsFunctor
{
   const Ray& ray_;
   real_t& t_;
   Vec3& n_;
   
   IntersectsFunctor(const Ray& ray, real_t& t, Vec3& n) : ray_(ray), t_(t), n_(n) {}
   
   template< typename BodyType >
   bool operator()( BodyType* bd1 ) {
      return intersects( bd1, ray_, t_, n_);
   }
};

inline bool intersects(const SphereID sphere, const Ray& ray, real_t& t, Vec3& n) {
   real_t inf = std::numeric_limits<real_t>::max();
   const Vec3& direction = ray.getDirection();
   Vec3 displacement = ray.getOrigin() - sphere->getPosition();
   real_t a = direction * direction;
   real_t b = real_t(2.) * (displacement * direction);
   real_t c = (displacement * displacement) - (sphere->getRadius() * sphere->getRadius());
   real_t discriminant = b*b - real_t(4.)*a*c;
   if (discriminant < EPSILON) {
      // with discriminant smaller than 0, sphere is not hit by ray
      // (no solution for quadratic equation)
      // with discriminant being 0, sphere only tangentially hits the ray (not enough)
      t = inf;
      return false;
   }
   real_t root = real_t(std::sqrt(discriminant));
   real_t t0 = (-b - root) / (real_t(2.) * a); // point where the ray enters the sphere
   real_t t1 = (-b + root) / (real_t(2.) * a); // point where the ray leaves the sphere
   if (t0 < 0 && t1 < 0) {
      t = inf;
      return false;
   }
   t = (t0 < t1) ? t0 : t1; // assign the closest hit point to t
   if (t < 0) {
      // at least one of the calculated hit points is behind the rays origin
      if (t1 < 0) {
         // both of the points are behind the origin (ray does not hit sphere)
         t = inf;
         return false;
      } else {
         // one point is hit by the ray (ray is within sphere)
         t = t1;
      }
   }
   Vec3 intersectionPoint = ray.getOrigin() + direction*t;
   n = (intersectionPoint - sphere->getPosition()).getNormalized();

   return true;
}

inline bool intersects(const PlaneID plane, const Ray& ray, real_t& t, Vec3& n) {
   real_t inf = std::numeric_limits<real_t>::max();
   const Vec3& direction = ray.getDirection();
   const Vec3& origin = ray.getOrigin();
   const Vec3& planeNormal = plane->getNormal();
   real_t denominator = planeNormal * direction;
   if (std::abs(denominator) > EPSILON) {
      real_t t_;
      t_ = ((plane->getPosition() - origin) * planeNormal) / denominator;
      if (t_ > EPSILON) {
         t = t_;
         n = planeNormal * sign(-denominator);
         return true;
      } else {
         t = inf;
      }
   }
   t = inf;
   return false;
}

inline bool intersects(const BoxID box, const Ray& ray, real_t& t, Vec3& n) {
   Ray transformedRay(box->pointFromWFtoBF(ray.getOrigin()), box->vectorFromWFtoBF(ray.getDirection()));
   
   const Vec3& lengths = box->getLengths();
   const Vec3 lengthsHalved = lengths/real_t(2);
   
   const Vec3 bounds[2] = {
      -lengthsHalved,
      lengthsHalved
   };
   
   const Vector3<int8_t>& sign = transformedRay.getInvDirectionSigns();
   const Vec3& invDirection = transformedRay.getInvDirection();
   const Vec3& origin = transformedRay.getOrigin();
   
   real_t inf = std::numeric_limits<real_t>::max();
   
   size_t tminAxis = 0, tmaxAxis = 0;
   real_t txmin, txmax;
   real_t tmin = txmin = (bounds[sign[0]][0] - origin[0]) * invDirection[0];
   real_t tmax = txmax = (bounds[1-sign[0]][0] - origin[0]) * invDirection[0];
   real_t tymin = (bounds[sign[1]][1] - origin[1]) * invDirection[1];
   real_t tymax = (bounds[1-sign[1]][1] - origin[1]) * invDirection[1];
   if (tmin > tymax || tymin > tmax) {
      t = inf;
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
      t = inf;
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
   
   n[0] = n[1] = n[2] = real_t(0);
   real_t t_;
   if (tmin > 0) {
      // ray hit box from outside
      t_ = tmin;
      n[tminAxis] = real_t(1);
   } else if (tmax < 0) {
      // tmin and tmax are smaller than 0 -> box is in rays negative direction
      t = inf;
      return false;
   } else {
      // ray origin within box
      t_ = tmax;
      n[tmaxAxis] = real_t(1);
   }
   
   if (transformedRay.getDirection() * n > 0) {
      n = -n;
   }
   
   n = box->vectorFromBFtoWF(n);
   WALBERLA_LOG_INFO("t_: " << t_ << ", n: " << n);
   
   t = t_;
   return true;
}

inline bool intersects(const AABB& aabb, const Ray& ray, real_t& t, real_t padding) {
   // An Efficient and Robust Ray–Box Intersection Algorithm: http://people.csail.mit.edu/amy/papers/box-jgt.pdf
   const Vec3 paddingVector(padding, padding, padding);
   Vec3 bounds[2] = {
      aabb.min() - paddingVector,
      aabb.max() + paddingVector
   };
   
   const Vector3<int8_t>& sign = ray.getInvDirectionSigns();
   const Vec3& invDirection = ray.getInvDirection();
   const Vec3& origin = ray.getOrigin();

   real_t inf = std::numeric_limits<real_t>::max();
   
   real_t txmin, txmax;
   real_t tmin = txmin = (bounds[sign[0]][0] - origin[0]) * invDirection[0];
   real_t tmax = txmax = (bounds[1-sign[0]][0] - origin[0]) * invDirection[0];
   real_t tymin = (bounds[sign[1]][1] - origin[1]) * invDirection[1];
   real_t tymax = (bounds[1-sign[1]][1] - origin[1]) * invDirection[1];
   if (tmin > tymax || tymin > tmax) {
      t = inf;
      return false;
   }
   if (tymin > tmin) {
      tmin = tymin;
   }
   if (tymax < tmax) {
      tmax = tymax;
   }
   real_t tzmin = (bounds[sign[2]][2] - origin[2]) * invDirection[2];
   real_t tzmax = (bounds[1-sign[2]][2] - origin[2]) * invDirection[2];
   if (tmin > tzmax || tzmin > tmax) {
      t = inf;
      return false;
   }
   if (tzmin > tmin) {
      tmin = tzmin;
   }
   if (tzmax < tmax) {
      tmax = tzmax;
   }
   
   real_t t_;
   if (tmin > 0) {
      // ray hit box from outside
      t_ = tmin;
   } else if (tmax < 0) {
      // tmin and tmax are smaller than 0 -> box is in rays negative direction
      t = inf;
      return false;
   } else {
      // ray origin within box
      t_ = tmax;
   }
   
   t = t_;
   return true;
}
}
}
}
