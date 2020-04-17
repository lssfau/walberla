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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesh_common/MatrixVectorOperations.h>
#include <mesh/pe/rigid_body/ConvexPolyhedron.h>
#include <mesh_common/TriangleMeshes.h>
#include <pe/raytracing/Intersects.h>

namespace walberla {
namespace pe {
namespace raytracing {

// Implemented following the description in
// "FAST RAY - CONVEX POLYHEDRON INTERSECTION" by Eric Haines
inline bool intersects(const mesh::pe::ConvexPolyhedronID poly, const Ray& ray, real_t& t_near, Vec3& n)
{
   Ray transformedRay = ray.transformedToBF(poly);

          t_near = std::numeric_limits<real_t>::min();
   real_t t_far  = std::numeric_limits<real_t>::max();
   const mesh::TriangleMesh& mesh = poly->getMesh();

   for(auto fh : mesh.faces())
   {
      //plane equation x*Pn + d = 0
      const Vector3<real_t>    Pn = mesh::toWalberla(mesh.normal(fh)); // Plane normal
      const real_t              d = - Pn * mesh::toWalberla(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(fh))));
      const real_t             vn = Pn * transformedRay.getOrigin() + d;
      const real_t             vd = Pn * transformedRay.getDirection();

      if ( floatIsEqual(vd, real_t(0)) )
      {
         if (vn > real_t(0)) return false;
         continue;
      }

      const real_t                 t = -vn / vd;

      if (vd > real_t(0))
      {
         // back-facing
         if (t < real_t(0)) return false;
         if (t < t_far) t_far=t;
      } else
      {
         // front-facing
         if (t > t_near)
         {
            t_near = t;
            n = Pn;
         }
         if (t_near > t_far) return false;
      }
   }

   if (t_near < 0) return false;

   n = poly->vectorFromBFtoWF(n);
   return true;
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
