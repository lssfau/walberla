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
//! \file Utility.cpp
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Helper functions for surface geometry computations.
//
//======================================================================================================================

#include "Utility.h"

#include "core/math/Constants.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"

#include <cmath>
#include <vector>

#include "ContactAngle.h"

namespace walberla
{
namespace free_surface
{
bool computeArtificalWallPoint(const Vector3< real_t >& globalInterfacePointLocation,
                               const Vector3< real_t >& globalCellCoordinate, const Vector3< real_t >& normal,
                               const Vector3< real_t >& wallNormal, const Vector3< real_t >& obstacleNormal,
                               const ContactAngle& contactAngle, Vector3< real_t >& artificialWallPointCoord,
                               Vector3< real_t >& artificialWallNormal)
{
   // get local interface point location (location of interface point inside cell with origin (0.5,0.5,0.5))
   const Vector3< real_t > interfacePointLocation = globalInterfacePointLocation - globalCellCoordinate;

   // check whether the interface plane intersects one of the cell's edges; exit this function if it does not intersect
   // any edge in at least one direction (see dissertation of S. Donath, 2011, section 6.4.5.4)
   if (wallNormal[0] < real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   if (wallNormal[0] > real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   if (wallNormal[1] < real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   if (wallNormal[1] > real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(1)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   if (wallNormal[2] < real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   if (wallNormal[2] > real_c(0) &&
       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(1)),
                                Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)) &&

       (getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(1)),
                                Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal,
                                interfacePointLocation) < real_c(0)))
   {
      return false;
   }

   // line 1 in Algorithm 6.2 in dissertation of S. Donath, 2011
   real_t cosAlpha = dot(normal, obstacleNormal);

   // determine sin(alpha) via orthogonal projections to compute alpha in the correct quadrant
   const Vector3< real_t > projector = normal - cosAlpha * obstacleNormal;
   const Vector3< real_t > baseVec   = cross(obstacleNormal, normal);
   const real_t sinAlpha             = projector * cross(baseVec, obstacleNormal).getNormalized();

   // compute alpha (angle between surface plane and wall)
   real_t alpha;
   if (sinAlpha >= real_c(0)) { alpha = real_c(std::acos(cosAlpha)); }
   else { alpha = real_c(2) * math::pi - real_c(std::acos(cosAlpha)); }

   // line 2 in Algorithm 6.2 in dissertation of S. Donath, 2011
   const real_t theta = contactAngle.getInRadians();
   const real_t delta = theta - alpha;

   // determine distance from interface point to wall plane
   const real_t wallDistance = dot(fabs(interfacePointLocation), wallNormal);

   // line 3-4 in Algorithm 6.2 in dissertation of S. Donath, 2011
   // correct contact angle is already reached
   if (realIsEqual(delta, real_c(0), real_c(1e-14)))
   {
      // IMPORTANT: the following approach is in contrast to the dissertation of S. Donath, 2011
      // - Dissertation: only the intersection of the interface surface with the wall is computed; it is known that this
      // could in rare cases lead to unstable configurations
      // - Here: extrapolate (extend) the wall point such that the triangle is considered valid during curvature
      // computation;
      // This has been adapted from the old waLBerla source code. While delta is considered to be zero, there is
      // still a division by delta below. This has not been found problematic when using double precision, however it
      // leads to invalid values in single precision. Therefore, the following macro exits the function prematurely when
      // using single precision and returns false such that the wall point will not be used.
#ifndef WALBERLA_DOUBLE_ACCURACY
      return false;
#endif

      // determine the direction in which the artifical wall point is to be expected
      // expression "dot(wallNormal,normal)*wallNormal" is orthogonal projection of normal on wallNormal
      // (wallNormal has already been normalized => no division by vector length required)
      const Vector3< real_t > targetDir = (normal - dot(wallNormal, normal) * wallNormal).getNormalized();

      const real_t sinTheta = contactAngle.getSin();

      // distance of interface point to wall in direction of wallNormal
      real_t wallPointDistance = wallDistance / sinTheta; // d_WP in dissertation of S. Donath, 2011

      // wall point must not be too close or far away from interface point too avoid degenerated triangles
      real_t virtWallDistance;
      if (wallPointDistance < real_c(0.8))
      {
         // extend distance with a heuristically chosen tolerance such that the point is not thrown away when checking
         // for degenerated triangles in the curvature computation
         wallPointDistance = real_c(0.801);
         virtWallDistance  = wallPointDistance * sinTheta; // d_W'
      }
      else
      {
         if (wallPointDistance > real_c(1.8))
         {
            // reduce distance with heuristically chosen tolerance
            wallPointDistance = real_c(1.799);
            virtWallDistance  = wallPointDistance * sinTheta; // d_W'
         }
         else
         {
            virtWallDistance = wallDistance; // d_W'
         }
      }

      const real_t virtPointDistance = virtWallDistance / sinTheta; // d_WP'

      // compute point by shifting virtWallDistance along wallNormal starting from globalInterfacePointLocation
      // virtWallProjection is given in global coordinates
      const Vector3< real_t > virtWallProjection = globalInterfacePointLocation - virtWallDistance * wallNormal;

      const real_t cosTheta = contactAngle.getCos();

      // compute artificial wall point by starting from virtWallProjection and shifting "virtPointDistance*cosTheta" in
      // targetDir (vritualWallPointCoord is r_W in dissertation of S. Donath, 2011)
      artificialWallPointCoord = virtWallProjection + virtPointDistance * cosTheta * targetDir;

      // radius r of the artificial circle "virtPointDistance*0.5/(sin(0.5)*delta)"
      // midpoint M of the artificial circle "normal*r+globalInterfacePointLocation"
      // normal of the virtual wall point "M-artificialWallPointCoord"
      artificialWallNormal = normal * virtPointDistance * real_c(0.5) / (std::sin(real_c(0.5) * delta)) +
                             globalInterfacePointLocation - artificialWallPointCoord;
      artificialWallNormal = artificialWallNormal.getNormalized();

      return true;
   }
   else
   {
      // compute base angles of triangle; line 6 in Algorithm 6.2 in dissertation of S. Donath, 2011
      const real_t beta = real_c(0.5) * (math::pi - real_c(std::fabs(delta)));

      real_t gamma;
      // line 7 in Algorithm 6.2 in dissertation of S. Donath, 2011
      if (theta < alpha) { gamma = beta - theta; }
      else { gamma = beta + theta - math::pi; }

      const real_t wallPointDistance =
         wallDistance / real_c(std::cos(std::fabs(gamma))); // d_WP in dissertation of S. Donath, 2011

      // line 9 in Algorithm 6.2 in dissertation of S. Donath, 2011
      // division by zero not possible as delta==0 is caught above
      real_t radius = real_c(0.5) * wallPointDistance / std::sin(real_c(0.5) * real_c(std::fabs(delta)));

      // check wallPointDistance for being in a valid range (to avoid degenerated triangles)
      real_t artificialWallPointDistance = wallPointDistance; // d'_WP in dissertation of S. Donath, 2011

      if (wallPointDistance < real_c(0.8))
      {
         // extend distance with a heuristically chosen tolerance such that the point is not thrown away when checking
         // for degenerated triangles in the curvature computation
         artificialWallPointDistance = real_c(0.801);

         // if extended distance exceeds circle diameter, assume delta=90 degrees
         if (artificialWallPointDistance > real_c(2) * radius)
         {
            radius = artificialWallPointDistance * math::one_div_root_two;
         }
      }
      else
      {
         // reduce distance with heuristically chosen tolerance
         if (wallPointDistance > real_c(1.8)) { artificialWallPointDistance = real_c(1.799); }
      }

      // line 17 in Algorithm 6.2 in dissertation of S. Donath, 2011
      const real_t artificialDelta =
         real_c(2) * real_c(std::asin(real_c(0.5) * artificialWallPointDistance /
                                      real_c(radius))); // delta' in dissertation of S. Donath, 2011

      // line 18 in Algorithm 6.2 in dissertation of S. Donath, 2011
      Vector3< real_t > rotVec = cross(normal, obstacleNormal);

      // change direction of rotation axis for delta>0; this is in contrast to Algorithm 6.2 in dissertation of Stefan
      // Donath, 2011 but was found to be necessary as otherwise the artificialWallNormal points in the wrong direction
      if (delta > real_c(0)) { rotVec *= real_c(-1); }

      // line 19 in Algorithm 6.2 in dissertation of S. Donath, 2011
      const Matrix3< real_t > rotMat(rotVec, artificialDelta);

      // line 20 in Algorithm 6.2 in dissertation of S. Donath, 2011
      artificialWallNormal = rotMat * normal;

      // line 21 in Algorithm 6.2 in dissertation of S. Donath, 2011
      if (theta < alpha)
      {
         artificialWallPointCoord = globalInterfacePointLocation + radius * (normal - artificialWallNormal);
      }
      else { artificialWallPointCoord = globalInterfacePointLocation - radius * (normal - artificialWallNormal); }

      return true;
   }
}

Vector3< real_t > getInterfacePoint(const Vector3< real_t >& normal, real_t fillLevel)
{
   // exploit symmetries of cubic cell to simplify the algorithm, i.e., restrict the fill level to the interval
   // [0,0.5] (see dissertation of T. Pohl, 2008, p. 22f)
   bool fillMirrored = false;
   if (fillLevel >= real_c(0.5))
   {
      fillMirrored = true;
      fillLevel    = real_c(1) - fillLevel;
   }

   // sort normal components such that nx, ny, nz >= 0 and nx <= ny <= nz to simplify the algorithm
   Vector3< real_t > normalSorted = fabs(normal);
   if (normalSorted[0] > normalSorted[1])
   {
      // swap nx and ny
      real_t tmp      = normalSorted[0];
      normalSorted[0] = normalSorted[1];
      normalSorted[1] = tmp;
   }

   if (normalSorted[1] > normalSorted[2])
   {
      // swap ny and nz
      real_t tmp      = normalSorted[1];
      normalSorted[1] = normalSorted[2];
      normalSorted[2] = tmp;

      if (normalSorted[0] > normalSorted[1])
      {
         // swap nx and ny
         tmp             = normalSorted[0];
         normalSorted[0] = normalSorted[1];
         normalSorted[1] = tmp;
      }
   }

   // minimal and maximal plane offset chosen as in the dissertation of T. Pohl, 2008, p. 22
   real_t offsetMin = real_c(-0.866025); // sqrt(3/4), lowest possible value
   real_t offsetMax = real_c(0);

   // find correct interface position by bisection (Algorithm 2.1, p. 22 in dissertation of T. Pohl, 2008)
   const uint_t numBisections =
      uint_c(10); // number of bisections, =10 in dissertation of T. Pohl, 2008 (heuristically chosen)

   for (uint_t i = uint_c(0); i <= numBisections; ++i)
   {
      const real_t offsetTmp    = real_c(0.5) * (offsetMin + offsetMax);
      const real_t newFillLevel = computeCellFluidVolume(normalSorted, offsetTmp);

      // volume is too small, reduce upper bound
      if (newFillLevel > fillLevel) { offsetMax = offsetTmp; }

      // volume is too large, reduce lower bound
      else { offsetMin = offsetTmp; }
   }
   real_t offset = real_c(0.5) * (offsetMin + offsetMax);

   if (fillMirrored) { offset *= real_c(-1); }
   const Vector3< real_t > interfacePoint = Vector3< real_t >(real_c(0.5)) + offset * normal;

   return interfacePoint;
}

real_t getCellEdgeIntersection(const Vector3< real_t >& edgePoint, const Vector3< real_t >& edgeDirection,
                               const Vector3< real_t >& normal, const Vector3< real_t >& surfacePoint)
{
   //#ifndef BELOW_CELL
   //#   define BELOW_CELL (-10)
   //#endif

#ifndef ABOVE_CELL
#   define ABOVE_CELL (-20)
#endif
   // mathematical description:
   // surface plane in coordinate form: x * normal = surfacePoint * normal
   // cell edge line: edgePoint + lambda * edgeDirection
   // => point of intersection: lambda = (interfacePoint * normal - edgePoint * normal) / edgeDirection * normal

   // compute angle between normal and cell edge
   real_t cosAngle = dot(edgeDirection, normal);

   real_t intersection = real_c(0);

   // intersection exists only if angle is not 90°, i.e., (intersection != 0)
   if (std::fabs(cosAngle) >= real_c(1e-14))
   {
      intersection = ((surfacePoint - edgePoint) * normal) / cosAngle;

      //      // intersection is below cell
      //      if (intersection < real_c(0)) { intersection = real_c(BELOW_CELL); }

      // intersection is above cell
      if (intersection > real_c(1)) { intersection = real_c(ABOVE_CELL); }
   }
   else // no intersection if angle is 90° (intersection == 0)
   {
      intersection = real_c(-1);
   }

   return intersection;
}

real_t computeCellFluidVolume(const Vector3< real_t >& normal, real_t offset)
{
   const Vector3< real_t > interfacePoint = Vector3< real_t >(real_c(0.5)) + offset * normal;

   real_t volume = real_c(0);

   // stores points of intersection with cell edges; points are shifted along normal such that the points lay
   // on one plane and surface area can be calculated
   std::vector< Vector3< real_t > > points;

   // SW: south west, EB: east bottom, etc.
   real_t iSW = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                        Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal, interfacePoint);

   // if no intersection with edge SW, volume is zero
   if (iSW > real_c(0) || realIsIdentical(iSW, ABOVE_CELL))
   {
      real_t iSE = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                           Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal, interfacePoint);
      real_t iNW = getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                           Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal, interfacePoint);
      real_t iNE = getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(1), real_c(0)),
                                           Vector3< real_t >(real_c(0), real_c(0), real_c(1)), normal, interfacePoint);

      // simple case: all four vertices are included in fluid domain (see Figure 2.12, p. 24 in dissertation of Thomas
      // Pohl, 2008)
      if (iNE > real_c(0)) { volume = real_c(0.25) * (iSW + iSE + iNW + iNE); }
      else
      {
         if (iSE >= real_c(0))
         {
            real_t iEB =
               getCellEdgeIntersection(Vector3< real_t >(real_c(1), real_c(0), real_c(0)),
                                       Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal, interfacePoint);

            // shift intersection points along normal and store points
            points.push_back(Vector3< real_t >(real_c(1), real_c(0), real_c(iSE)) - interfacePoint);
            points.push_back(Vector3< real_t >(real_c(1), real_c(iEB), real_c(0)) - interfacePoint);

            volume += iSE * iEB;
         }
         else
         {
            real_t iSB =
               getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                       Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal, interfacePoint);

            points.push_back(Vector3< real_t >(real_c(iSB), real_c(0), real_c(0)) - interfacePoint);
         }

         if (iNW >= real_c(0))
         {
            real_t iNB =
               getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(1), real_c(0)),
                                       Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal, interfacePoint);

            points.push_back(Vector3< real_t >(real_c(iNB), real_c(1), real_c(0)) - interfacePoint);
            points.push_back(Vector3< real_t >(real_c(0), real_c(1), real_c(iNW)) - interfacePoint);

            volume += iNB * iNW;
         }
         else
         {
            real_t iWB =
               getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(0)),
                                       Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal, interfacePoint);

            points.push_back(Vector3< real_t >(real_c(0), real_c(iWB), real_c(0)) - interfacePoint);
         }

         real_t iWT =
            getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                    Vector3< real_t >(real_c(0), real_c(1), real_c(0)), normal, interfacePoint);
         if (iWT >= real_c(0))
         {
            real_t iST =
               getCellEdgeIntersection(Vector3< real_t >(real_c(0), real_c(0), real_c(1)),
                                       Vector3< real_t >(real_c(1), real_c(0), real_c(0)), normal, interfacePoint);

            points.push_back(Vector3< real_t >(real_c(0), real_c(iWT), real_c(1)) - interfacePoint);
            points.push_back(Vector3< real_t >(real_c(iST), real_c(0), real_c(1)) - interfacePoint);

            volume += iWT * iST;
         }
         else { points.push_back(Vector3< real_t >(real_c(0), real_c(0), real_c(iSW)) - interfacePoint); }

         Vector3< real_t > vectorProduct(real_c(0));
         real_t area = real_c(0);
         size_t j    = points.size() - 1;

         // compute the vector product, the length of which gives the surface area of the corresponding
         // parallelogram;
         for (size_t i = 0; i != points.size(); ++i)
         {
            vectorProduct[0] = points[i][1] * points[j][2] - points[i][2] * points[j][1];
            vectorProduct[1] = points[i][2] * points[j][0] - points[i][0] * points[j][2];
            vectorProduct[2] = points[i][0] * points[j][1] - points[i][1] * points[j][0];

            // area of the triangle is obtained through division by 2; this is done below in the volume
            // calculation (where the division is then by 6 instead of by 3)
            area += vectorProduct.length();
            j = i;
         }

         // compute and sum the volumes of each resulting pyramid: V = area * height * 1/3
         volume += area * (normal * interfacePoint);
         volume /= real_c(6.0); // division by 6 since the area was not divided by 2 above
      }
   }
   else // no intersection with edge SW, volume is zero
   {
      volume = real_c(0);
   }

   return volume;
}
} // namespace free_surface
} // namespace walberla
