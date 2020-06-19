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
//! \file ContactDetection.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/collision_detection/AnalyticCollisionFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/BaseShape.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

class ContactDetection : public collision_detection::AnalyticContactDetection
{
public:
   ///import all collision detection functions from AnalyticContactDetection
   using AnalyticContactDetection::operator();

   ///overwrite functions that should be different
   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::Sphere& geo1,
                    const data::Sphere& geo2,
                    Accessor& ac);

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::Sphere& s,
                    const data::HalfSpace& p,
                    Accessor& ac );
};

template <typename Accessor>
inline bool ContactDetection::operator()( const size_t idx1,
                                          const size_t idx2,
                                          const data::Sphere& geo1,
                                          const data::Sphere& geo2,
                                          Accessor& ac)
{
   using namespace collision_detection::analytic;
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   //ensure collision order idx2 has to be larger than idx1
   if (ac.getUid(idx2) < ac.getUid(idx1))
      return operator()(idx2, idx1, geo2, geo1, ac);

   getIdx1() = idx1;
   getIdx2() = idx2;
   return detectSphereSphereCollision(ac.getPosition(getIdx1()),
                                      ac.getRadiusAtTemperature(getIdx1()),
                                      ac.getPosition(getIdx2()),
                                      ac.getRadiusAtTemperature(getIdx2()),
                                      getContactPoint(),
                                      getContactNormal(),
                                      getPenetrationDepth(),
                                      getContactThreshold());
}

template <typename Accessor>
inline bool ContactDetection::operator()( const size_t idx1,
                                          const size_t idx2,
                                          const data::Sphere& s,
                                          const data::HalfSpace& p,
                                          Accessor& ac )
{
   using namespace collision_detection::analytic;
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   getIdx1() = idx1;
   getIdx2() = idx2;
   return detectSphereHalfSpaceCollision(ac.getPosition(getIdx1()),
                                         ac.getRadiusAtTemperature(getIdx1()),
                                         ac.getPosition(getIdx2()),
                                         p.getNormal(),
                                         getContactPoint(),
                                         getContactNormal(),
                                         getPenetrationDepth(),
                                         getContactThreshold());
}

} //namespace mesa_pd
} //namespace walberla
