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
#include "pe/rigidbody/Ellipsoid.h"
#include "pe/rigidbody/Union.h"
#include "pe/utility/BodyCast.h"
#include <core/math/Utility.h>

#include <pe/raytracing/Ray.h>

namespace walberla {
namespace pe {
namespace raytracing {

inline bool intersects(const SphereID sphere, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const PlaneID plane, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const BoxID box, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const CapsuleID capsule, const Ray& ray, real_t& t, Vec3& n);
inline bool intersects(const EllipsoidID ellipsoid, const Ray& ray, real_t& t, Vec3& n);

inline bool intersects(const BodyID body, const Ray& ray, real_t& t, Vec3& n);
   
inline bool intersects(const AABB& aabb, const Ray& ray, real_t& t, real_t padding = real_t(0.0), Vec3* n = NULL);
inline bool intersectsSphere(const Vec3& gpos, real_t radius, const Ray& ray, real_t& t0, real_t& t1);

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

const real_t discriminantEps = real_t(1e-4);
   
inline bool intersects(const SphereID sphere, const Ray& ray, real_t& t, Vec3& n) {
   const real_t realMax = std::numeric_limits<real_t>::max();
   
   real_t t0, t1;
   if (!intersectsSphere(sphere->getPosition(), sphere->getRadius(), ray, t0, t1)) {
      t = realMax;
      return false;
   }
   
   t = (t0 < t1) ? t0 : t1; // assign the closest hit point to t
   if (t < 0) {
      t = t1; // if t0 < 0, t1 is > 0 (else intersectsSphere would have returned false)
   }
   
   Vec3 intersectionPoint = ray.getOrigin() + ray.getDirection()*t;
   n = (intersectionPoint - sphere->getPosition()).getNormalized();

   return true;
}

inline bool intersects(const PlaneID plane, const Ray& ray, real_t& t, Vec3& n) {
   const real_t realMax = std::numeric_limits<real_t>::max();
   const Vec3& direction = ray.getDirection();
   const Vec3& origin = ray.getOrigin();
   const Vec3& planeNormal = plane->getNormal();
   real_t denominator = planeNormal * direction;
   if (std::abs(denominator) > discriminantEps) {
      real_t t_;
      t_ = ((plane->getPosition() - origin) * planeNormal) / denominator;
      if (t_ > 0) {
         t = t_;
         n = planeNormal * walberla::math::sign(-denominator);
         return true;
      } else {
         t = realMax;
      }
   }
   t = realMax;
   return false;
}

inline bool intersects(const BoxID box, const Ray& ray, real_t& t, Vec3& n) {
   Ray transformedRay = ray.transformedToBF(box);
   
   const Vec3& lengths = box->getLengths();
   const Vec3 lengthsHalved = lengths/real_t(2);
   
   const Vec3 bounds[2] = {
      -lengthsHalved,
      lengthsHalved
   };
   
   const Vector3<int8_t>& sign = transformedRay.getInvDirectionSigns();
   const Vec3& invDirection = transformedRay.getInvDirection();
   const Vec3& origin = transformedRay.getOrigin();
   
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
   
   n[0] = n[1] = n[2] = real_t(0);
   real_t t_;
   if (tmin > 0) {
      // ray hit box from outside
      t_ = tmin;
      n[tminAxis] = real_t(1);
   } else if (tmax < 0) {
      // tmin and tmax are smaller than 0 -> box is in rays negative direction
      t = realMax;
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
   
   t = t_;
   return true;
}
   
inline bool intersects(const CapsuleID capsule, const Ray& ray, real_t& t, Vec3& n) {
   const Ray transformedRay = ray.transformedToBF(capsule);
   const Vec3& direction = transformedRay.getDirection();
   const Vec3& origin = transformedRay.getOrigin();
   real_t halfLength = capsule->getLength()/real_t(2);

   const real_t realMax = std::numeric_limits<real_t>::max();
   t = realMax;
   
   bool t0hit = false, t1hit = false;
   size_t intersectedPrimitive = 0; // 1 for capsule, 2 for left half sphere, 3 for right half sphere

   real_t a = direction[2]*direction[2] + direction[1]*direction[1];
   real_t b = real_t(2)*origin[2]*direction[2] + real_t(2)*origin[1]*direction[1];
   real_t c = origin[2]*origin[2] + origin[1]*origin[1] - capsule->getRadius()*capsule->getRadius();
   real_t discriminant = b*b - real_t(4.)*a*c;
   if (std::abs(discriminant) >= discriminantEps) {
      // With discriminant smaller than 0, cylinder is not hit by ray (no solution for quadratic equation).
      // Thus only enter this section if the equation is actually solvable.
      
      real_t root = real_t(std::sqrt(discriminant));
      real_t t0 = (-b - root) / (real_t(2) * a); // Distance to point where the ray enters the cylinder
      real_t t1 = (-b + root) / (real_t(2) * a); // Distance to point where the ray leaves the cylinder

      real_t tx0 = origin[0] + direction[0]*t0;
      real_t tx1 = origin[0] + direction[0]*t1;
      
      if (t0 > 0 && tx0 >= -halfLength && tx0 <= halfLength) {
         t0hit = true;
         intersectedPrimitive = 1;
         t = t0;
      }
      if (t1 > 0 && tx1 >= -halfLength && tx1 <= halfLength && t1 < t) {
         t1hit = true;
         if (t1 < t) {
            intersectedPrimitive = 1;
            t = t1;
         }
      }
   }
   
   // Check now for end capping half spheres.
   // Only check them if the ray didnt both enter and leave the cylinder part of the capsule already (t0hit && t1hit).
   if (!t0hit || !t1hit) {
      real_t t0_left, t1_left;
      Vec3 leftSpherePos(-halfLength, 0, 0);
      if (intersectsSphere(leftSpherePos, capsule->getRadius(), transformedRay, t0_left, t1_left)) {
         // at least one of t0_left and t1_left are not behind the rays origin
         real_t t0x_left = origin[0] + direction[0]*t0_left;
         real_t t1x_left = origin[0] + direction[0]*t1_left;
         
         real_t t_left = realMax;
         if (t0_left > 0 && t0x_left < -halfLength) {
            t_left = t0_left;
         }
         if (t1_left > 0 && t1x_left < -halfLength && t1_left < t_left) {
            t_left = t1_left;
         }
         if (t_left < t) {
            intersectedPrimitive = 2;
            t = t_left;
         }
      }
      
      real_t t0_right, t1_right;
      Vec3 rightSpherePos(halfLength, 0, 0);
      if (intersectsSphere(rightSpherePos, capsule->getRadius(), transformedRay, t0_right, t1_right)) {
         // At least one of t0_right and t1_right are not behind the rays origin
         real_t t0x_right = origin[0] + direction[0]*t0_right;
         real_t t1x_right = origin[0] + direction[0]*t1_right;
         
         real_t t_right = realMax;
         if (t0_right > 0 && t0x_right > halfLength) {
            t_right = t0_right;
         }
         if (t1_right > 0 && t1x_right > halfLength && t1_right < t_right) {
            t_right = t1_right;
         }
         if (t_right < t) {
            intersectedPrimitive = 3;
            t = t_right;
         }
      }
      
      if (realIsIdentical(t, realMax)) {
         return false;
      }
      
      Vec3 intersectionPoint = origin + direction*t;
      if (intersectedPrimitive == 2) {
         n = (intersectionPoint - leftSpherePos).getNormalized();
      } else if (intersectedPrimitive == 3) {
         n = (intersectionPoint - rightSpherePos).getNormalized();
      }
   }
   
   WALBERLA_ASSERT(intersectedPrimitive != 0);
   
   if (intersectedPrimitive == 1) {
      Vec3 intersectionPoint = origin + direction*t;
      Vec3 intersectionPointOnXAxis(intersectionPoint[0], real_t(0), real_t(0));
      n = (intersectionPoint - intersectionPointOnXAxis).getNormalized();
   }
   
   n = capsule->vectorFromBFtoWF(n);
   
   return true;
}
   
inline bool intersects(const EllipsoidID ellipsoid, const Ray& ray, real_t& t, Vec3& n) {
   const real_t realMax = std::numeric_limits<real_t>::max();

   const Ray transformedRay = ray.transformedToBF(ellipsoid);
   const Vec3& semiAxes = ellipsoid->getSemiAxes();
   
   const Mat3 M = Mat3::makeDiagonalMatrix(real_t(1)/semiAxes[0], real_t(1)/semiAxes[1], real_t(1)/semiAxes[2]);
   
   const Vec3 d_M = M*transformedRay.getDirection();
   const Vec3 P_M = M*transformedRay.getOrigin();
   
   const real_t a = d_M*d_M;
   const real_t b = real_t(2)*P_M*d_M;
   const real_t c = P_M*P_M - 1;
   
   const real_t discriminant = b*b - real_t(4.)*a*c;
   if (discriminant < 0) {
      // with discriminant smaller than 0, sphere is not hit by ray
      // (no solution for quadratic equation)
      t = realMax;
      return false;
   }
   
   const real_t root = real_t(std::sqrt(discriminant));
   const real_t t0 = (-b - root) / (real_t(2.) * a); // distance to point where the ray enters the sphere
   const real_t t1 = (-b + root) / (real_t(2.) * a); // distance to point where the ray leaves the sphere
   
   if (t0 < 0 && t1 < 0) {
      return false;
   }
   t = (t0 < t1) ? t0 : t1; // assign the closest distance to t
   if (t < 0) {
      // at least one of the calculated distances is behind the rays origin
      if (t1 < 0) {
         // both of the points are behind the origin (ray does not hit sphere)
         return false;
      } else {
         t = t1;
      }
   }
   
   const Vec3 transformedN = transformedRay.getOrigin() + t*transformedRay.getDirection();
   const Mat3 M_inv = M.getInverse();
   n = ellipsoid->vectorFromBFtoWF((M_inv*transformedN).getNormalized());
   
   return true;
}
   
inline bool intersects(const BodyID body, const Ray& ray, real_t& t, Vec3& n) {
   WALBERLA_UNUSED(body);
   WALBERLA_UNUSED(ray);
   WALBERLA_UNUSED(t);
   WALBERLA_UNUSED(n);
   WALBERLA_ABORT("This ray - body intersection test is not implemented yet!");
   return false;
}
   
inline bool intersectsSphere(const Vec3& gpos, real_t radius, const Ray& ray, real_t& t0, real_t& t1) {
   const real_t realMax = std::numeric_limits<real_t>::max();
   
   const Vec3& direction = ray.getDirection();
   Vec3 displacement = ray.getOrigin() - gpos;
   
   real_t a = direction * direction;
   real_t b = real_t(2.) * (displacement * direction);
   real_t c = (displacement * displacement) - (radius * radius);
   real_t discriminant = b*b - real_t(4.)*a*c;
   if (discriminant < 0) {
      // with discriminant smaller than 0, sphere is not hit by ray
      // (no solution for quadratic equation)
      t0 = realMax;
      t1 = realMax;
      return false;
   }
   
   real_t root = real_t(std::sqrt(discriminant));
   t0 = (-b - root) / (real_t(2.) * a); // distance to point where the ray enters the sphere
   t1 = (-b + root) / (real_t(2.) * a); // distance to point where the ray leaves the sphere
   
   if (t0 < 0 && t1 < 0) {
      return false;
   }
   real_t t = (t0 < t1) ? t0 : t1; // assign the closest distance to t
   if (t < 0) {
      // at least one of the calculated distances is behind the rays origin
      if (t1 < 0) {
         // both of the points are behind the origin (ray does not hit sphere)
         return false;
      }
   }
   
   return true;
}

inline bool intersects(const AABB& aabb, const Ray& ray, real_t& t, real_t padding, Vec3* n) {
   // An Efficient and Robust Ray–Box Intersection Algorithm: http://people.csail.mit.edu/amy/papers/box-jgt.pdf
   const Vec3 paddingVector(padding, padding, padding);
   Vec3 bounds[2] = {
      aabb.min() - paddingVector,
      aabb.max() + paddingVector
   };
   
   const Vector3<int8_t>& sign = ray.getInvDirectionSigns();
   const Vec3& invDirection = ray.getInvDirection();
   const Vec3& origin = ray.getOrigin();

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
   
   if (n != NULL) {
      (*n)[0] = (*n)[1] = (*n)[2] = real_t(0);
   }
   real_t t_;
   if (tmin > 0) {
      // ray hit box from outside
      t_ = tmin;
      if (n != NULL) {
         (*n)[tminAxis] = real_t(1);
      }
   } else if (tmax < 0) {
      // tmin and tmax are smaller than 0 -> box is in rays negative direction
      t = realMax;
      return false;
   } else {
      // ray origin within box
      t_ = tmax;
      if (n != NULL) {
         (*n)[tmaxAxis] = real_t(1);
      }
   }
   
   if (n != NULL) {
      if (ray.getDirection() * (*n) > 0) {
         *n = -(*n);
      }
   }
   
   t = t_;
   return true;
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
