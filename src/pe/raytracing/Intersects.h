#pragma once

#include <pe/Types.h>
#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/CylindricalBoundary.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"
#include "pe/utility/BodyCast.h"

#include <pe/raytracing/Ray.h>

#include <boost/tuple/tuple.hpp>

#define EPSILON real_t(1e-4)

namespace walberla {
   namespace pe {
      namespace raytracing {
         inline bool intersects(SphereID sphere, Ray* ray, real_t &t);
         inline bool intersects(PlaneID plane, Ray* ray, real_t* t);
         inline bool intersects(BoxID box, Ray* ray, real_t* t);
         inline bool intersects(GeomPrimitive* body, Ray* ray, real_t* t);
         
         struct IntersectsFunctor
         {
            Ray* ray_;
            real_t* t_;
            
            IntersectsFunctor(Ray* ray, real_t* t) : ray_(ray), t_(t) {}
            
            template< typename BodyType >
            bool operator()( BodyType* bd1 ) { return intersects( bd1, ray_, t_); }
         };
         
         inline bool intersects(SphereID sphere, Ray* ray, real_t* t) {
            Vec3 displacement = ray->origin - sphere->getPosition();
            real_t a = ray->direction * ray->direction;
            real_t b = real_t(2.) * (displacement * ray->direction);
            real_t c = (displacement * displacement) - (sphere->getRadius() * sphere->getRadius());
            real_t discriminant = b*b - real_t(4.)*a*c;
            if (discriminant < EPSILON) {
               // with discriminant smaller than 0, sphere is not hit by ray
               // (no solution for quadratic equation)
               // with discriminant being 0, sphere only tangentially hits the ray (not enough)
               *t = real_t(INFINITY);
               return false;
            }
            real_t root = real_t(std::sqrt(discriminant));
            real_t t0 = (-b - root) / (real_t(2.) * a); // point where the ray enters the sphere
            real_t t1 = (-b + root) / (real_t(2.) * a); // point where the ray leaves the sphere
            if (t0 < 0 && t1 < 0) {
               *t = real_t(INFINITY);
               return false;
            }
            real_t t_;
            t_ = (t0 < t1) ? t0 : t1; // assign the closest hit point to t
            if (t_ < 0) {
               // at least one of the calculated hit points is behind the rays origin
               if (t1 < 0) {
                  // both of the points are behind the origin (ray does not hit sphere)
                  *t = real_t(INFINITY);
                  return false;
               } else {
                  // one point is hit by the ray (ray is within sphere)
                  t_ = t1;
               }
            }
            *t = t_;
            return true;
         }
         
         inline bool intersects(PlaneID plane, Ray* ray, real_t* t) {
            real_t denominator = plane->getNormal() * ray->direction;
            if (std::abs(denominator) > EPSILON) {
               real_t t_;
               t_ = ((plane->getPosition() - ray->origin) * plane->getNormal()) / denominator;
               
               if (t_ > EPSILON) {
                  Vec3 intersectionPoint = (t_*ray->direction + ray->origin);
                  Vec3 originToIntersection = ray->origin - intersectionPoint;
                  *t = originToIntersection.length();
                  return true;
               } else {
                  *t = INFINITY;
               }
            }
            *t = INFINITY;
            return false;
         }
         
         inline bool intersects(BoxID box, Ray* ray, real_t* t) {
            Mat3 invRotationMatrix = box->getRotation().getInverse();
            Ray transformedRay(invRotationMatrix * (ray->origin - box->getPosition()), invRotationMatrix * ray->direction);
            
            Vec3 lengths = box->getLengths();
            Vec3 lengthsHalved = lengths/real_t(2);
            
            Vec3 bounds[2] = {
               -lengthsHalved,
               lengthsHalved
            };
                        
            real_t txmin, txmax;
            real_t tmin = txmin = (bounds[transformedRay.sign[0]][0] - transformedRay.origin[0]) * transformedRay.inv_direction[0];
            real_t tmax = txmax = (bounds[1-transformedRay.sign[0]][0] - transformedRay.origin[0]) * transformedRay.inv_direction[0];
            real_t tymin = (bounds[transformedRay.sign[1]][1] - transformedRay.origin[1]) * transformedRay.inv_direction[1];
            real_t tymax = (bounds[1-transformedRay.sign[1]][1] - transformedRay.origin[1]) * transformedRay.inv_direction[1];
            if (tmin > tymax || tymin > tmax) {
               *t = INFINITY;
               return false;
            }
            if (tymin > tmin) {
               tmin = tymin;
            }
            if (tymax < tmax) {
               tmax = tymax;
            }
            real_t tzmin = (bounds[transformedRay.sign[2]][2] - transformedRay.origin[2]) * transformedRay.inv_direction[2];
            real_t tzmax = (bounds[1-transformedRay.sign[2]][2] - transformedRay.origin[2]) * transformedRay.inv_direction[2];
            if (tmin > tzmax || tzmin > tmax) {
               *t = INFINITY;
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
               *t = INFINITY;
               return false;
            } else {
               // ray origin within box
               t_ = tmax;
            }
            
            *t = t_;
            return true;
         }
         
         inline bool intersects(GeomPrimitive* body, Ray* ray, real_t* t) {
            // An Efficient and Robust Ray–Box Intersection Algorithm: http://people.csail.mit.edu/amy/papers/box-jgt.pdf
            AABB aabb = body->getAABB();
            
            Vec3 bounds[2] = {
               aabb.min(),
               aabb.max()
            };
            
            real_t txmin, txmax;
            real_t tmin = txmin = (bounds[ray->sign[0]][0] - ray->origin[0]) * ray->inv_direction[0];
            real_t tmax = txmax = (bounds[1-ray->sign[0]][0] - ray->origin[0]) * ray->inv_direction[0];
            real_t tymin = (bounds[ray->sign[1]][1] - ray->origin[1]) * ray->inv_direction[1];
            real_t tymax = (bounds[1-ray->sign[1]][1] - ray->origin[1]) * ray->inv_direction[1];
            if (tmin > tymax || tymin > tmax) {
               *t = INFINITY;
               return false;
            }
            if (tymin > tmin) {
               tmin = tymin;
            }
            if (tymax < tmax) {
               tmax = tymax;
            }
            real_t tzmin = (bounds[ray->sign[2]][2] - ray->origin[2]) * ray->inv_direction[2];
            real_t tzmax = (bounds[1-ray->sign[2]][2] - ray->origin[2]) * ray->inv_direction[2];
            if (tmin > tzmax || tzmin > tmax) {
               *t = INFINITY;
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
               *t = INFINITY;
               return false;
            } else {
               // ray origin within box
               t_ = tmax;
            }
            
            *t = t_;
            return true;
         }
      }
   }
}
