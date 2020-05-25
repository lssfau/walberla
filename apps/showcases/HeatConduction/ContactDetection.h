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

class ContactDetection
{
public:
   size_t& getIdx1() {return idx1_;}
   size_t& getIdx2() {return idx2_;}
   Vec3&   getContactPoint() {return contactPoint_;}
   Vec3&   getContactNormal() {return contactNormal_;}
   real_t& getPenetrationDepth() {return penetrationDepth_;}

   const size_t& getIdx1() const {return idx1_;}
   const size_t& getIdx2() const {return idx2_;}
   const Vec3&   getContactPoint() const {return contactPoint_;}
   const Vec3&   getContactNormal() const {return contactNormal_;}
   const real_t& getPenetrationDepth() const {return penetrationDepth_;}

   real_t& getContactThreshold() {return contactThreshold_;}

   template <typename GEO1_T, typename GEO2_T, typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const GEO1_T& geo1,
                    const GEO2_T& geo2,
                    Accessor& ac);

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

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::HalfSpace& p,
                    const data::Sphere& s,
                    Accessor& ac);

private:
   size_t idx1_;
   size_t idx2_;
   Vec3   contactPoint_;
   Vec3   contactNormal_;
   real_t penetrationDepth_;

   real_t contactThreshold_ = real_t(0.0);
};

template <typename GEO1_T, typename GEO2_T, typename Accessor>
inline bool ContactDetection::operator()( const size_t /*idx1*/,
                                          const size_t /*idx2*/,
                                          const GEO1_T& /*geo1*/,
                                          const GEO2_T& /*geo2*/,
                                          Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision not implemented!")
}

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
                                      ac.getRadius(getIdx1()),
                                      ac.getPosition(getIdx2()),
                                      ac.getRadius(getIdx2()),
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
                                         ac.getRadius(getIdx1()),
                                         ac.getPosition(getIdx2()),
                                         p.getNormal(),
                                         getContactPoint(),
                                         getContactNormal(),
                                         getPenetrationDepth(),
                                         getContactThreshold());
}

template <typename Accessor>
inline bool ContactDetection::operator()( const size_t idx1,
                                          const size_t idx2,
                                          const data::HalfSpace& p,
                                          const data::Sphere& s,
                                          Accessor& ac)
{
   return operator()(idx2, idx1, s, p, ac);
}

inline
std::ostream& operator<<( std::ostream& os, const ContactDetection& ac )
{
   os << "idx1:               " << ac.getIdx1() << "\n" <<
         "idx2:               " << ac.getIdx2() << "\n" <<
         "contact point:      " << ac.getContactPoint() << "\n" <<
         "contact normal:     " << ac.getContactNormal() << "\n" <<
         "penetration depth:  " << ac.getPenetrationDepth() << std::endl;
   return os;
}

} //namespace mesa_pd
} //namespace walberla
