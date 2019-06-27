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
//! \file GeneralContactDetection.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/collision_detection/EPA.h>
#include <mesa_pd/collision_detection/GJK.h>
#include <mesa_pd/collision_detection/Support.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/BaseShape.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

/**
 * Collision detection functor which uses recommended collision functions.
 *
 * Calculates and stores contact information between two particles.
 * If a collision was successfully detected by the operator()
 * contactPoint, contactNormal and penetrationDepth contain the
 * contact information. If no collision was detected the values of
 * these variables is undefined!
 */
class GeneralContactDetection : public AnalyticContactDetection
{
public:
   using AnalyticContactDetection::operator();

   template <typename GEO1_T, typename GEO2_T, typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const GEO1_T& geo1,
                    const GEO2_T& geo2,
                    Accessor& ac);

   template <typename GEO2_T, typename Accessor>
   bool operator()(const size_t idx1,
                   const size_t idx2,
                   const data::HalfSpace& geo1,
                   const GEO2_T& geo2,
                   Accessor& ac);

   template <typename GEO1_T, typename Accessor>
   bool operator()(const size_t idx1,
                   const size_t idx2,
                   const GEO1_T& geo1,
                   const data::HalfSpace& geo2,
                   Accessor& ac);

   template <typename GEO2_T, typename Accessor>
   bool operator()(const size_t idx1,
                   const size_t idx2,
                   const data::CylindricalBoundary& geo1,
                   const GEO2_T& geo2,
                   Accessor& ac);

   template <typename GEO1_T, typename Accessor>
   bool operator()(const size_t idx1,
                   const size_t idx2,
                   const GEO1_T& geo1,
                   const data::CylindricalBoundary& geo2,
                   Accessor& ac);
private:
   bool collideGJKEPA(Support& geom0, Support& geom1);
};

template <typename GEO1_T, typename GEO2_T, typename Accessor>
bool GeneralContactDetection::operator()( const size_t idx1,
                                          const size_t idx2,
                                          const GEO1_T& geo1,
                                          const GEO2_T& geo2,
                                          Accessor& ac)
{
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   //ensure collision order idx2 has to be larger than idx1
   if (ac.getUid(idx2) < ac.getUid(idx1))
      return operator()(idx2, idx1, geo2, geo1, ac);

   getIdx1() = idx1;
   getIdx2() = idx2;

   Support e1(ac.getPosition(getIdx1()), ac.getRotation(getIdx1()), geo1);
   Support e2(ac.getPosition(getIdx2()), ac.getRotation(getIdx2()), geo2);
   return collideGJKEPA(e1, e2);
}

template <typename GEO2_T, typename Accessor>
bool GeneralContactDetection::operator()(const size_t idx1,
                                         const size_t idx2,
                                         const data::HalfSpace& geo1,
                                         const GEO2_T& geo2,
                                         Accessor& ac)
{
   getIdx1() = idx1;
   getIdx2() = idx2;

   Support sup(ac.getPosition(idx2), ac.getRotation(idx2), geo2);
   Vec3 support_dir = -geo1.getNormal();
   // We now have a direction facing to the "wall".
   // Compute support point of body b in this direction. This will be the furthest point overlapping.
   Vec3 contactp = sup.support(support_dir);
   real_t pdepth = contactp * geo1.getNormal() - ac.getPosition(idx1) * geo1.getNormal();
   if(pdepth < contactThreshold_)
   { //We have a collision
      contactNormal_ = support_dir;
      penetrationDepth_ = pdepth;
      contactPoint_ = contactp + real_t(0.5) * penetrationDepth_ * contactNormal_;
      return true;
   } else
   { //No collision
      return false;
   }
}

template <typename GEO1_T, typename Accessor>
bool GeneralContactDetection::operator()(const size_t idx1,
                                         const size_t idx2,
                                         const GEO1_T& geo1,
                                         const data::HalfSpace& geo2,
                                         Accessor& ac)
{
   return operator()(idx2, idx1, geo2, geo1, ac);
}

template <typename GEO2_T, typename Accessor>
bool GeneralContactDetection::operator()(const size_t idx1,
                                         const size_t idx2,
                                         const data::CylindricalBoundary& geo1,
                                         const GEO2_T& geo2,
                                         Accessor& ac)
{
   getIdx1() = idx1;
   getIdx2() = idx2;

   Support sup(ac.getPosition(idx2), ac.getRotation(idx2), geo2);

   WALBERLA_CHECK_FLOAT_EQUAL(geo1.getAxis().sqrLength(), real_t(1));

   auto d = ac.getPosition(idx2) - (dot((ac.getPosition(idx2) - ac.getPosition(idx1)), geo1.getAxis()) * geo1.getAxis() + ac.getPosition(idx1));
   Vec3 farestP = sup.support(d);
   auto d2 = farestP - (dot((farestP - ac.getPosition(idx1)), geo1.getAxis()) * geo1.getAxis() + ac.getPosition(idx1));
   real_t dist = d2.sqrLength();

   if(dist > geo1.getRadius() * geo1.getRadius())
   { //We have a collision
      penetrationDepth_ = geo1.getRadius() - std::sqrt(dist);
      normalize(d2);
      contactNormal_ = d2;
      contactPoint_ = farestP + d2 * penetrationDepth_ * real_t(0.5);
      return true;
   } else
   { //No collision
      return false;
   }
}

template <typename GEO1_T, typename Accessor>
bool GeneralContactDetection::operator()(const size_t idx1,
                                         const size_t idx2,
                                         const GEO1_T& geo1,
                                         const data::CylindricalBoundary& geo2,
                                         Accessor& ac)
{
   return operator()(idx2, idx1, geo2, geo1, ac);
}

bool GeneralContactDetection::collideGJKEPA(Support& geom0, Support& geom1)
{
   real_t margin = real_t(1e-4);
   GJK gjk;
   if(gjk.doGJKmargin(geom0, geom1, margin))
   {
      EPA epa;
      epa.useSphereOptimization(false);
      if (epa.doEPAmargin(geom0, geom1, gjk, contactNormal_, contactPoint_, penetrationDepth_, margin))
      {
         return true;
      } else
      {
         return false;
      }
   } else
   {
      return false;
   }
}

} //namespace collision_detection
} //namespace mesa_pd
} //namespace walberla
