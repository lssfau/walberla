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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/collision_detection/AnalyticCollisionFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/BaseShape.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

/**
 * Collision detection functor which uses analytic functions.
 *
 * Calculates and stores contact information between two particles.
 * If a collision was successfully detected by the operator()
 * contactPoint, contactNormal and penetrationDepth contain the
 * contact information. If no collision was detected the values of
 * these variables is undefined!
 */
class AnalyticContactDetection
{
public:
   virtual ~AnalyticContactDetection() = default;

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
                    const data::Box& bx,
                    Accessor& ac );

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::Box& bx,
                    const data::Sphere& s,
                    Accessor& ac);

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::Sphere& s,
                    const data::CylindricalBoundary& cb,
                    Accessor& ac );

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::CylindricalBoundary& cb,
                    const data::Sphere& s,
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

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::HalfSpace& p0,
                    const data::HalfSpace& p1,
                    Accessor& ac);

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::CylindricalBoundary& p0,
                    const data::HalfSpace& p1,
                    Accessor& ac);

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::HalfSpace& p0,
                    const data::CylindricalBoundary& p1,
                    Accessor& ac);

   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    const data::CylindricalBoundary& p0,
                    const data::CylindricalBoundary& p1,
                    Accessor& ac);

protected:
   size_t idx1_ = std::numeric_limits<size_t>::max();
   size_t idx2_ = std::numeric_limits<size_t>::max();
   Vec3   contactPoint_ = Vec3();
   Vec3   contactNormal_ = Vec3();
   real_t penetrationDepth_ = real_t(0);

   real_t contactThreshold_ = real_t(0.0);
};

template <typename GEO1_T, typename GEO2_T, typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t /*idx1*/,
                                                  const size_t /*idx2*/,
                                                  const GEO1_T& /*geo1*/,
                                                  const GEO2_T& /*geo2*/,
                                                  Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision not implemented!")
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::Sphere& geo1,
                                                  const data::Sphere& geo2,
                                                  Accessor& ac)
{
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   //ensure collision order idx2 has to be larger than idx1
   if (ac.getUid(idx2) < ac.getUid(idx1))
      return operator()(idx2, idx1, geo2, geo1, ac);

   getIdx1() = idx1;
   getIdx2() = idx2;
   return analytic::detectSphereSphereCollision(ac.getPosition(getIdx1()),
                                                geo1.getRadius(),
                                                ac.getPosition(getIdx2()),
                                                geo2.getRadius(),
                                                getContactPoint(),
                                                getContactNormal(),
                                                getPenetrationDepth(),
                                                getContactThreshold());
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::Sphere& s,
                                                  const data::Box& bx,
                                                  Accessor& ac )
{
   getIdx1() = idx1;
   getIdx2() = idx2;
   return analytic::detectSphereBoxCollision(ac.getPosition(getIdx1()),
                                             s.getRadius(),
                                             ac.getPosition(getIdx2()),
                                             bx.getEdgeLength(),
                                             ac.getRotation(getIdx2()),
                                             getContactPoint(),
                                             getContactNormal(),
                                             getPenetrationDepth(),
                                             getContactThreshold());
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::Box& bx,
                                                  const data::Sphere& s,
                                                  Accessor& ac)
{
   return operator()(idx2, idx1, s, bx, ac);
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::Sphere& s,
                                                  const data::CylindricalBoundary& cb,
                                                  Accessor& ac )
{
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   getIdx1() = idx1;
   getIdx2() = idx2;
   return analytic::detectSphereCylindricalBoundaryCollision(ac.getPosition(getIdx1()),
                                                             s.getRadius(),
                                                             ac.getPosition(getIdx2()),
                                                             cb.getRadius(),
                                                             cb.getAxis(),
                                                             getContactPoint(),
                                                             getContactNormal(),
                                                             getPenetrationDepth(),
                                                             getContactThreshold());
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::CylindricalBoundary& cb,
                                                  const data::Sphere& s,
                                                  Accessor& ac)
{
   return operator()(idx2, idx1, s, cb, ac);
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::Sphere& s,
                                                  const data::HalfSpace& p,
                                                  Accessor& ac )
{
   WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "colliding with itself!");

   getIdx1() = idx1;
   getIdx2() = idx2;
   return analytic::detectSphereHalfSpaceCollision(ac.getPosition(getIdx1()),
                                                   s.getRadius(),
                                                   ac.getPosition(getIdx2()),
                                                   p.getNormal(),
                                                   getContactPoint(),
                                                   getContactNormal(),
                                                   getPenetrationDepth(),
                                                   getContactThreshold());
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t idx1,
                                                  const size_t idx2,
                                                  const data::HalfSpace& p,
                                                  const data::Sphere& s,
                                                  Accessor& ac)
{
   return operator()(idx2, idx1, s, p, ac);
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t /*idx1*/,
                                                  const size_t /*idx2*/,
                                                  const data::HalfSpace& /*geo1*/,
                                                  const data::HalfSpace& /*geo2*/,
                                                  Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision between two half spaces is not defined!")
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t /*idx1*/,
                                                  const size_t /*idx2*/,
                                                  const data::HalfSpace& /*geo1*/,
                                                  const data::CylindricalBoundary& /*geo2*/,
                                                  Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision between half spaces and cylindrical boundary is not defined!")
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t /*idx1*/,
                                                  const size_t /*idx2*/,
                                                  const data::CylindricalBoundary& /*geo1*/,
                                                  const data::HalfSpace& /*geo2*/,
                                                  Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision between half spaces and cylindrical boundary is not defined!")
}

template <typename Accessor>
inline bool AnalyticContactDetection::operator()( const size_t /*idx1*/,
                                                  const size_t /*idx2*/,
                                                  const data::CylindricalBoundary& /*geo1*/,
                                                  const data::CylindricalBoundary& /*geo2*/,
                                                  Accessor& /*ac*/)
{
   WALBERLA_ABORT("Collision between two cylindrical boundaries is not defined!")
}

inline
std::ostream& operator<<( std::ostream& os, const AnalyticContactDetection& ac )
{
   os << "idx1:               " << ac.getIdx1() << "\n" <<
         "idx2:               " << ac.getIdx2() << "\n" <<
         "contact point:      " << ac.getContactPoint() << "\n" <<
         "contact normal:     " << ac.getContactNormal() << "\n" <<
         "penetration depth:  " << ac.getPenetrationDepth() << std::endl;
   return os;
}

} //namespace collision_detection
} //namespace mesa_pd
} //namespace walberla
