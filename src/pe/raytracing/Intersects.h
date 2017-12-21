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
         inline bool intersects(PlaneID plane, Ray* ray, real_t &t);
         
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
            real_t root = sqrt(discriminant);
            real_t t0 = (-b - root) / (real_t(2.) * a); // point where the ray enters the sphere
            real_t t1 = (-b + root) / (real_t(2.) * a); // point where the ray leaves the sphere
            if (t0 < 0 && t1 < 0) {
               *t = real_t(INFINITY);
               return false;
            }
            real_t t_;
            t_ = (t0 < t1) ? t0 : t1; // assign the closest hit point to t
            if (t_ < 0) {
               // at least one of the calculated hit point is behind the rays origin
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
      }
   }
}
