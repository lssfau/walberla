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
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/IAccessor.h>

#include <mesa_pd/common/ParticleFunctions.h>

#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include <mesa_pd/data/shape/ConvexPolyhedron.h>

namespace walberla {
namespace mesa_pd {

/*
 * "contains point" functionality
 * can either be formulated in world frame coordinates (then the rotation of the geometry is not taken into account)
 * or in body frame coordinates (BF) which requires the point to be first transformed
 */

inline
bool isPointInsideSphere(const Vec3& point,
                         const Vec3& spherePosition, const real_t sphereRadius )
{
   return !((point - spherePosition).sqrLength() > sphereRadius * sphereRadius);
}

inline
bool isPointInsideHalfSpace(const Vec3& point,
                            const Vec3& halfSpacePosition, const Vec3& halfSpaceNormal )
{
   return !((point - halfSpacePosition) * halfSpaceNormal > real_t(0));
}

inline
bool isPointInsideBoxBF(const Vec3& pointBF,
                        const Vec3& edgeLengths )
{
   return std::fabs(pointBF[0]) <= real_t(0.5)*edgeLengths[0] &&
          std::fabs(pointBF[1]) <= real_t(0.5)*edgeLengths[1] &&
          std::fabs(pointBF[2]) <= real_t(0.5)*edgeLengths[2];
}

inline
bool isPointInsideEllipsoidBF(const Vec3& pointBF,
                              const Vec3& semiAxes )
{
   return ( (pointBF[0] * pointBF[0])/(semiAxes[0] * semiAxes[0]) + (pointBF[1] * pointBF[1])/(semiAxes[1] * semiAxes[1])
            + (pointBF[2] * pointBF[2])/(semiAxes[2] * semiAxes[2]) <= 1_r );
}

inline
bool isPointInsideCylindricalBoundary(const Vec3& point,
                                      const Vec3& cylindricalBoundaryPosition, const real_t radius, const Vec3& axis  )
{
   const Vec3 distanceFromCylinderCenterLine = (point - cylindricalBoundaryPosition) - ((point - cylindricalBoundaryPosition) * axis) * axis;
   return distanceFromCylinderCenterLine.sqrLength() >= radius*radius;
}

#ifdef WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE
inline
bool isPointInsideConvexPolyhedronBF(const Vec3& point, const mesh::TriangleMesh& mesh)
{
   WALBERLA_ASSERT(mesh.has_face_normals(), "Provided mesh has no face normals! E.g., call `mesh.request_face_normals(); mesh.update_face_normals();` to add them.")

   return std::none_of(mesh.faces().begin(),
                       mesh.faces().end(),
                       [&](auto fh)
                       {
                         //check if point is on positive side of the face
                         const mesh::TriangleMesh::Normal &n = mesh.normal(fh); // Plane normal
                         const mesh::TriangleMesh::Point &pp = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(fh))); // Point on plane

                         // normal * (p - pointOnPlane)
                         return (n[0] * (point[0] - pp[0]) + n[1] * (point[1] - pp[1]) + n[2] * (point[2] - pp[2]) >= real_t(0));
                       });
}
#endif

struct ContainsPointFunctor
{
   template<typename ParticleAccessor_T, typename Shape_T>
   bool operator()(const size_t /*particleIdx*/, const Shape_T& /*shape*/, const ParticleAccessor_T& /*ac*/, const Vector3<real_t>& /*point*/ )
   {
      WALBERLA_ABORT("ContainsPoint not implemented!");
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::Sphere& sphere, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideSphere(point, ac.getPosition(particleIdx), sphere.getRadius() );
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::HalfSpace& halfSpace, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideHalfSpace(point, ac.getPosition(particleIdx), halfSpace.getNormal() );
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::Box& box, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideBoxBF(mesa_pd::transformPositionFromWFtoBF(particleIdx, ac, point), box.getEdgeLength());
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::Ellipsoid& ellipsoid, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideEllipsoidBF(mesa_pd::transformPositionFromWFtoBF(particleIdx, ac, point), ellipsoid.getSemiAxes());
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::CylindricalBoundary& cylindricalBoundary, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideCylindricalBoundary(point, ac.getPosition(particleIdx), cylindricalBoundary.getRadius(), cylindricalBoundary.getAxis());
   }

#ifdef WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE
   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::ConvexPolyhedron& convexPolyhedron, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      auto pointBF = mesa_pd::transformPositionFromWFtoBF(particleIdx, ac, point);
      if( pointBF.sqrLength() > convexPolyhedron.getBoundingSphereRadius() * convexPolyhedron.getBoundingSphereRadius() )
         return false;

      return isPointInsideConvexPolyhedronBF(pointBF, convexPolyhedron.getMesh());
   }
#endif

};


} //namespace mesa_pd
} //namespace walberla
