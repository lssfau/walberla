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
//! \file CurvatureSweep.impl.h
//! \ingroup surface_geometry
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweeps for computing the interface curvature.
//
//======================================================================================================================

#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Matrix3.h"
#include "core/math/Utility.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include "stencil/D3Q27.h"
#include "stencil/Directions.h"

#include <algorithm>
#include <cmath>

#include "ContactAngle.h"
#include "CurvatureSweep.h"
#include "Utility.h"

namespace walberla
{
namespace free_surface
{
template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void CurvatureSweepFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(
   IBlock* const block)
{
   // get fields
   ScalarField_T* const curvatureField            = block->getData< ScalarField_T >(curvatureFieldID_);
   const VectorField_T* const normalField         = block->getData< const VectorField_T >(normalFieldID_);
   const VectorField_T* const obstacleNormalField = block->getData< const VectorField_T >(obstacleNormalFieldID_);
   const FlagField_T* const flagField             = block->getData< const FlagField_T >(flagFieldID_);

   // get flags
   const flag_t interfaceFlag              = flagField->getFlag(interfaceFlagID_);
   const flag_t liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);
   const flag_t obstacleFlagMask           = flagField->getMask(obstacleFlagIDSet_);

   WALBERLA_FOR_ALL_CELLS(
      flagFieldIt, flagField, normalFieldIt, normalField, obstacleNormalFieldIt, obstacleNormalField, curvatureFieldIt,
      curvatureField, {
         real_t& curv = *curvatureFieldIt;
         curv         = real_c(0.0);

         if (isFlagSet(flagFieldIt, interfaceFlag)) // only treat interface cells
         {
            real_t weightSum = real_c(0);

            if (normalFieldIt->sqrLength() < real_c(1e-14))
            {
               WALBERLA_LOG_WARNING("Invalid normal detected in CurvatureSweep.")
               continue;
            }

            // Parker-Youngs central finite difference approximation of curvature (see dissertation
            // of S. Bogner, 2017, section 4.4.2.1)
            for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
            {
               const Vector3< real_t > dirVector =
                  Vector3< real_t >(real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()));

               Vector3< real_t > neighborNormal;

               if (isPartOfMaskSet(flagFieldIt.neighbor(*dir), liquidInterfaceGasFlagMask | obstacleFlagMask))
               {
                  // get interface normal of neighbor in direction dir with respect to wetting model
                  if (enableWetting_ && isPartOfMaskSet(flagFieldIt.neighbor(*dir), obstacleFlagMask))
                  {
                     neighborNormal = getNormalWithWetting(normalFieldIt, obstacleNormalFieldIt, *dir);
                     neighborNormal = neighborNormal.getNormalizedOrZero();
                  }
                  else
                  {
                     if (isPartOfMaskSet(flagFieldIt.neighbor(*dir), liquidInterfaceGasFlagMask))
                     {
                        neighborNormal = normalFieldIt.neighbor(*dir);
                     }
                     else
                     {
                        // skip remainder of this direction such that it is not considered in curvature computation
                        continue;
                     }
                  }
               }

               // equation (35) in Brackbill et al. discretized with finite difference method
               // according to Parker-Youngs
               if constexpr (Stencil_T::D == uint_t(2))
               {
                  const real_t weight =
                     real_c(stencil::gaussianMultipliers[stencil::D3Q27::idx[stencil::map2Dto3D[2][*dir]]]);
                  weightSum += weight;
                  curv += weight * (dirVector * neighborNormal);
               }
               else
               {
                  const real_t weight = real_c(stencil::gaussianMultipliers[dir.toIdx()]);
                  weightSum += weight;
                  curv += weight * (dirVector * neighborNormal);
               }
            }

            // divide by sum of weights in Parker-Youngs approximation; sum does not contain weights of directions in
            // which there is no liquid, interface, gas, or obstacle cell (only when wetting is enabled, otherwise
            // obstacle cell is also not considered); must be done like this because otherwise such non-valid
            // directions would implicitly influence the finite difference scheme by assuming a normal of zero
            curv /= weightSum;
         }
      }) // WALBERLA_FOR_ALL_CELLS
}

template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void CurvatureSweepLocalTriangulation< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(
   IBlock* const block)
{
   const auto blockForest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   // struct for storing relevant information of each neighboring interface cell (POD type)
   using Neighbor = struct
   {
      Vector3< real_t > diff;     // difference (distance in coordinates) between this and neighboring interface point
      Vector3< real_t > diffNorm; // normalized difference between this and neighboring interface point
      Vector3< real_t > normal;   // interface normal of neighboring interface point
      real_t dist2;               // square of the distance between this and neighboring interface point
      real_t area;                // sum of the area of the two triangles that this neighboring interface is part of
      real_t sort;                // angle that is used to sort the order neighboring points accordingly
      bool valid;                 // validity, used to remove triangles with too narrow angles
      bool wall;                  // used to include wetting effects near solid cells
   };

   // get fields
   ScalarField_T* const curvatureField            = block->getData< ScalarField_T >(curvatureFieldID_);
   const VectorField_T* const normalField         = block->getData< const VectorField_T >(normalFieldID_);
   const ScalarField_T* const fillField           = block->getData< const ScalarField_T >(fillFieldID_);
   const FlagField_T* const flagField             = block->getData< const FlagField_T >(flagFieldID_);
   const VectorField_T* const obstacleNormalField = block->getData< const VectorField_T >(obstacleNormalFieldID_);

   // get flags
   auto interfaceFlag    = flagField->getFlag(interfaceFlagID_);
   auto obstacleFlagMask = flagField->getMask(obstacleFlagIDSet_);

   WALBERLA_FOR_ALL_CELLS(
      flagFieldIt, flagField, normalFieldIt, normalField, fillFieldIt, fillField, obstacleNormalFieldIt,
      obstacleNormalField, curvatureFieldIt, curvatureField, {
         real_t& curv = *curvatureFieldIt;
         curv         = real_c(0.0);
         std::vector< Neighbor > neighbors;
         Vector3< real_t > meanInterfaceNormal;

         // compute curvature only in interface points
         if (!isFlagSet(flagFieldIt, interfaceFlag)) { continue; }

         // normal of this cell also contributes to mean normal
         meanInterfaceNormal = *normalFieldIt;

         // iterate over all neighboring cells (Eq. 2.18 in dissertation of T. Pohl, 2008)
         for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
         {
            auto neighborFlags = flagFieldIt.neighbor(*dir);

            if (isFlagSet(neighborFlags, interfaceFlag) || // Eq. 2.19 in dissertation of T. Pohl, 2008
                (isPartOfMaskSet(neighborFlags, obstacleFlagMask) &&
                 dir.toIdx() <= uint_c(18))) // obstacle in main direction or diagonal direction (not corner direction)
            {
               Vector3< real_t > neighborNormal;
               const real_t neighborFillLevel = fillFieldIt.neighbor(*dir);

               // if loop was entered because of neighboring solid cell, normal of this solid cell points towards the
               // currently processed interface cell
               Vector3< real_t > wallNormal(real_c(-dir.cx()), real_c(-dir.cy()), real_c(-dir.cz()));

               if (isPartOfMaskSet(neighborFlags, obstacleFlagMask))
               {
                  neighborNormal =
                     Vector3< real_t >(real_c(1)); // temporarily guarantees "neighborNormal.sqrLength()>0"
                  wallNormal.getNormalized();
               }
               else { neighborNormal = normalFieldIt.neighbor(*dir); }

               if (neighborNormal.sqrLength() > real_c(0))
               {
                  Neighbor n;

                  // get global coordinate (with respect to whole simulation domain) of the currently processed cell
                  Vector3< real_t > globalCellCoordinate =
                     blockForest->getBlockLocalCellCenter(*block, flagFieldIt.cell()) - Vector3< real_t >(real_c(0.5));

                  for (auto dir2 = Stencil_T::beginNoCenter(); dir2 != Stencil_T::end(); ++dir2)
                  {
                     // stay in the close neighborhood of the currently processed interface cell
                     if ((dir.cx() != 0 && dir2.cx() != 0) || (dir.cy() != 0 && dir2.cy() != 0) ||
                         (dir.cz() != 0 && dir2.cz() != 0))
                     {
                        continue;
                     }

                     if (isPartOfMaskSet(neighborFlags, obstacleFlagMask) && enableWetting_)
                     {
                        // get flags of neighboring cell in direction dir2
                        auto neighborFlagsInDir2 = flagFieldIt.neighbor(*dir2);

                        // the currently processed interface cell i has a neighboring solid cell j in direction dir;
                        // get the flags of j's neighboring cell in direction dir2
                        // i.e., from the current cell, go to neighbor in dir; from there, go to next cell in dir2
                        auto neighborNeighborFlags =
                           flagFieldIt.neighbor(dir.cx() + dir2.cx(), dir.cy() + dir2.cy(), dir.cz() + dir2.cz());

                        // true if the currently processed interface cell i has a neighboring interface cell j (in
                        // dir2), which (j) has a neighboring obstacle cell in the same direction as i does (in dir)
                        if (isFlagSet(neighborFlagsInDir2, interfaceFlag) &&
                            isPartOfMaskSet(neighborNeighborFlags, obstacleFlagMask))
                        {
                           // get the normal of the currently processed interface cell's neighboring interface cell
                           // (in direction 2)
                           Vector3< real_t > neighborNormalDir2 = normalFieldIt.neighbor(*dir2);

                           Vector3< real_t > neighborGlobalCoordDir2 = globalCellCoordinate;
                           neighborGlobalCoordDir2[0] += real_c(dir2.cx());
                           neighborGlobalCoordDir2[1] += real_c(dir2.cy());
                           neighborGlobalCoordDir2[2] += real_c(dir2.cz());

                           // get neighboring interface point, i.e., location of interface within cell
                           Vector3< real_t > neighborGlobalInterfacePoint =
                              getInterfacePoint(normalFieldIt.neighbor(*dir2), fillFieldIt.neighbor(*dir2));

                           // transform to global coordinates, i.e., neighborInterfacePoint specifies the global
                           // location of the interface point in the currently processed interface cell's neighbor in
                           // direction dir2
                           neighborGlobalInterfacePoint += neighborGlobalCoordDir2;

                           // get the mean (averaged over multiple solid cells) wall normal of the neighbor in
                           // direction dir2
                           Vector3< real_t > obstacleNormal = obstacleNormalFieldIt.neighbor(*dir2);
                           obstacleNormal *= real_c(-1);
                           if (obstacleNormal.sqrLength() < real_c(1e-10)) { obstacleNormal = wallNormal; }

                           Vector3< real_t > neighborPoint;

                           bool result = computeArtificalWallPoint(
                              neighborGlobalInterfacePoint, neighborGlobalCoordDir2, neighborNormalDir2, wallNormal,
                              obstacleNormal, contactAngle_, neighborPoint, neighborNormal);
                           if (!result) { continue; }
                           n.wall  = true;
                           n.diff  = neighborPoint - neighborGlobalInterfacePoint;
                           n.dist2 = n.diff.sqrLength();
                        }
                        else { continue; }
                     }
                     else
                     // will be entered if:
                     // isFlagSet(neighborFlags, interfaceFlag) && !isPartOfMaskSet(neighborFlags, obstacleFlagMask)
                     {
                        n.wall = false;

                        // get neighboring interface point, i.e., location of interface within cell
                        n.diff = getInterfacePoint(neighborNormal, neighborFillLevel);

                        // get distance between this cell (0,0,0) and neighboring interface point + (dx,dy,dz)
                        n.diff += Vector3< real_t >(real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()));

                        // get distance between this and neighboring interface point
                        n.diff -= getInterfacePoint(*normalFieldIt, *fillFieldIt);
                        n.dist2 = n.diff.sqrLength();
                     }

                     // exclude neighboring interface points that are too close or too far away from this cell's
                     // interface point
                     if (n.dist2 >= real_c(0.64) && n.dist2 <= real_c(3.24)) // Eq. 2.20, 0.64 = 0.8^2; 3.24 = 1.8^2
                     {
                        n.normal   = neighborNormal;
                        n.diffNorm = n.diff.getNormalized();
                        n.area     = real_c(0);
                        n.sort     = real_c(0);
                        n.valid    = true;

                        neighbors.push_back(n);
                     }

                     // if there is no obstacle, loop should be interrupted immediately
                     if (!isPartOfMaskSet(neighborFlags, obstacleFlagMask))
                     {
                        // interrupt loop
                        break;
                     }
                  }
               }
            }
         }

         // remove degenerated triangles, see dissertation of T. Pohl, 2008, p. 27
         for (auto nIt1 = ++neighbors.begin(); !neighbors.empty() && nIt1 != neighbors.end(); ++nIt1)
         {
            if (!nIt1->valid) { continue; }

            for (auto nIt2 = neighbors.begin(); nIt2 != nIt1; ++nIt2)
            {
               if (!nIt2->valid) { continue; }

               // triangle is degenerated if angle between surface normals is less than 30° (heuristically chosen
               // in dissertation of T. Pohl, 2008, p. 27); greater sign is correct here due to cos(29) > cos(30)
               if (nIt1->diffNorm * nIt2->diffNorm >
                   real_c(0.866)) // cos(30°) = 0.866, as in dissertation of T. Pohl, 2008, p. 27
               {
                  const real_t diff = nIt1->dist2 - nIt2->dist2;

                  if (diff < real_c(1e-4)) { nIt1->valid = nIt1->wall ? true : false; }
                  if (diff > real_c(-1e-4)) { nIt2->valid = nIt2->wall ? true : false; }
               }
            }
         }

         // remove invalid neighbors
         neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(), [](const Neighbor& a) { return !a.valid; }),
                         neighbors.end());

         if (neighbors.size() < 4)
         {
            // WALBERLA_LOG_WARNING_ON_ROOT(
            //    "Not enough faces in curvature reconstruction, setting curvature in this cell to zero.");
            curv = real_c(0); // not documented in literature but taken from S. Donath's code
            continue;         // process next cell in WALBERLA_FOR_ALL_CELLS
         }

         // compute mean normal
         for (auto const& n : neighbors)
         {
            meanInterfaceNormal += n.normal;
         }
         meanInterfaceNormal = meanInterfaceNormal.getNormalized();

         // compute xAxis and yAxis that define a coordinate system on a tangent plane for sorting neighbors
         // T_i' = (I - N * N^t) * (p_i - p); projection of neighbors.diff[0] onto tangent plane (Figure 2.14 in
         // dissertation of T. Pohl, 2008)
         Vector3< real_t > xAxis =
            (Matrix3< real_t >::makeIdentityMatrix() - dyadicProduct(meanInterfaceNormal, meanInterfaceNormal)) *
            neighbors[0].diff;

         // T_i = T_i' / ||T_i'||
         xAxis = xAxis.getNormalized();

         // get vector that is orthogonal to xAxis and meanInterfaceNormal
         const Vector3< real_t > yAxis = cross(xAxis, meanInterfaceNormal);

         for (auto& n : neighbors)
         {
            // get cosine of angles between n.diff and axes of the new coordinate system
            const real_t cosAngX = xAxis * n.diff;
            const real_t cosAngY = yAxis * n.diff;

            // sort the neighboring interface points using atan2 (which is not just atan(wy/wx), see Wikipedia)
            n.sort = std::atan2(cosAngY, cosAngX);
         }

         std::sort(neighbors.begin(), neighbors.end(),
                   [](const Neighbor& a, const Neighbor& b) { return a.sort < b.sort; });

         Vector3< real_t > meanTriangleNormal(real_c(0));
         for (auto nIt1 = neighbors.begin(); neighbors.size() > uint_c(1) && nIt1 != neighbors.end(); ++nIt1)
         {
            // index of second neighbor starts over at 0: (k + 1) mod N_P
            auto nIt2 = (nIt1 != (neighbors.end() - 1)) ? (nIt1 + 1) : neighbors.begin();

            // N_f_k (with real length, i.e., not normalized yet)
            const Vector3< real_t > triangleNormal = cross(nIt1->diff, nIt2->diff);

            // |f_k| (multiplication with 0.5, since triangleNormal.length() is area of parallelogram and not
            // triangle)
            const real_t area = real_c(0.5) * triangleNormal.length();

            // lambda_i' = |f_{(i-1+N_p) mod N_p}| + |f_i|
            nIt1->area += area; // from the view of na, this is |f_i|
            nIt2->area += area; // from the view of nb, this is |f_{(i-1+N_p) mod N_p}|, i.e., area of the face
                                // from the neighbor with smaller index

            // N' = sum(|f_k| * N_f_k)
            meanTriangleNormal += area * triangleNormal;
         }

         if (meanTriangleNormal.length() < real_c(1e-10))
         {
            // WALBERLA_LOG_WARNING_ON_ROOT("Invalid meanTriangleNormal, setting curvature in this cell to zero.");
            curv = real_c(0); // not documented in literature but taken from S. Donath's code
            continue;         // process next cell in WALBERLA_FOR_ALL_CELLS
         }

         // N = N' / ||N'||
         meanTriangleNormal = meanTriangleNormal.getNormalized();

         // I - N * N^t; matrix for projection of vector on tangent plane
         const Matrix3< real_t > projMatrix =
            Matrix3< real_t >::makeIdentityMatrix() - dyadicProduct(meanTriangleNormal, meanTriangleNormal);

         // M
         Matrix3< real_t > mMatrix(real_c(0));

         // sum(lambda_i')
         real_t neighborAreaSum = real_c(0);

         // M = sum(lambda_i' * kappa_i * T_i * T_i^t)
         for (auto& n : neighbors)
         {
            if (n.area > real_c(0))
            {
               // kappa_i = 2 * N^t * (p_i - p) / ||p_i - p||^2
               const real_t kappa = real_c(2) * (meanTriangleNormal * n.diff) / n.dist2;

               // T_i' = (I - N * N^t) * (p_i - p)
               Vector3< real_t > tVector = projMatrix * n.diff;

               // T_i = T_i' / ||T_i'||
               tVector = tVector.getNormalized();

               // T_i * T_i^t
               const Matrix3< real_t > auxMat = dyadicProduct(tVector, tVector);

               // M += T_i * T_i^t * kappa_i * lambda_i'
               mMatrix += auxMat * kappa * n.area;

               // sum(lambda_i')
               neighborAreaSum += n.area;
            }
         }

         // M = M * 1 / sum(lambda_i')
         mMatrix = mMatrix * (real_c(1) / neighborAreaSum);

         // kappa = tr(M)
         curv = (mMatrix(0, 0) + mMatrix(1, 1) + mMatrix(2, 2));
      }) // WALBERLA_FOR_ALL_CELLS
}

template< typename Stencil_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void CurvatureSweepSimpleFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(
   IBlock* const block)
{
   // get fields
   ScalarField_T* const curvatureField    = block->getData< ScalarField_T >(curvatureFieldID_);
   const VectorField_T* const normalField = block->getData< const VectorField_T >(normalFieldID_);
   const FlagField_T* const flagField     = block->getData< const FlagField_T >(flagFieldID_);

   // get flags
   auto interfaceFlag    = flagField->getFlag(interfaceFlagID_);
   auto obstacleFlagMask = flagField->getMask(obstacleFlagIDSet_);

   WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, curvatureFieldIt, curvatureField, {
      real_t& curv = *curvatureFieldIt;
      curv         = real_c(0.0);

      // interface cells
      if (isFlagSet(flagFieldIt, interfaceFlag))
      {
         // this cell is next to a wall/obstacle
         if (enableWetting_ && isFlagInNeighborhood< Stencil_T >(flagFieldIt, obstacleFlagMask))
         {
            // compute wall/obstacle curvature
            vector_t obstacleNormal(real_c(0), real_c(0), real_c(0));
            for (auto it = Stencil_T::beginNoCenter(); it != Stencil_T::end(); ++it)
            {
               // calculate obstacle normal with central finite difference approximation of the surface's gradient (see
               // dissertation of S. Donath, 2011, section 6.3.5.2)
               if (isPartOfMaskSet(flagFieldIt.neighbor(*it), obstacleFlagMask))
               {
                  obstacleNormal[0] -= real_c(it.cx());
                  obstacleNormal[1] -= real_c(it.cy());
                  obstacleNormal[2] -= real_c(it.cz());
               }
            }

            if (obstacleNormal.sqrLength() > real_c(0))
            {
               obstacleNormal = obstacleNormal.getNormalized();

               // IMPORTANT REMARK:
               // the following wetting model is not documented in literature and not tested for correctness; use it
               // with caution
               curv = -real_c(0.25) * (contactAngle_.getCos() - (*normalFieldIt) * obstacleNormal);
            }
         }
         else // no obstacle cell is in next neighborhood
         {
            // central finite difference approximation of curvature (see dissertation of S. Bogner, 2017,
            // section 4.4.2.1)
            curv = normalFieldIt.neighbor(1, 0, 0)[0] - normalFieldIt.neighbor(-1, 0, 0)[0] +
                   normalFieldIt.neighbor(0, 1, 0)[1] - normalFieldIt.neighbor(0, -1, 0)[1] +
                   normalFieldIt.neighbor(0, 0, 1)[2] - normalFieldIt.neighbor(0, 0, -1)[2];

            curv *= real_c(0.25);
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

} // namespace free_surface
} // namespace walberla
