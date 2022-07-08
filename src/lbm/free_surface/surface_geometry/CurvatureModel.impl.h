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
//! \file CurvatureModel.impl.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Collection of sweeps required for using a specific curvature model.
//
//======================================================================================================================

#include "lbm/blockforest/communication/UpdateSecondGhostLayer.h"

#include "CurvatureModel.h"
#include "CurvatureSweep.h"
#include "DetectWettingSweep.h"
#include "ExtrapolateNormalsSweep.h"
#include "NormalSweep.h"
#include "ObstacleFillLevelSweep.h"
#include "ObstacleNormalSweep.h"
#include "SmoothingSweep.h"
#include "SurfaceGeometryHandler.h"

namespace walberla
{
namespace free_surface
{
namespace curvature_model
{
// empty sweep required for using selectors (e.g. StateSweep::fullFreeSurface)
struct emptySweep
{
   void operator()(IBlock*) {}
};

template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
void FiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::addSweeps(
   SweepTimeloop& timeloop, const FiniteDifferenceMethod::SurfaceGeometryHandler_T& geometryHandler)
{
   using Communication_T = typename SurfaceGeometryHandler_T::Communication_T;
   using StateSweep      = typename SurfaceGeometryHandler_T::StateSweep;

   // layout for allocating the smoothed fill level field
   field::Layout fillFieldLayout = field::fzyx;

   // check if an obstacle cell is in a non-periodic outermost global ghost layer; used to check if two ghost layers are
   // required for the fill level field
   const Vector3< bool > isObstacleInGlobalGhostLayerXYZ =
      geometryHandler.freeSurfaceBoundaryHandling_->isObstacleInGlobalGhostLayer();

   bool isObstacleInGlobalGhostLayer = false;
   if ((!geometryHandler.blockForest_->isXPeriodic() && isObstacleInGlobalGhostLayerXYZ[0]) ||
       (!geometryHandler.blockForest_->isYPeriodic() && isObstacleInGlobalGhostLayerXYZ[1]) ||
       (!geometryHandler.blockForest_->isZPeriodic() && isObstacleInGlobalGhostLayerXYZ[2]))
   {
      isObstacleInGlobalGhostLayer = true;
   }

   for (auto blockIt = geometryHandler.blockForest_->begin(); blockIt != geometryHandler.blockForest_->end(); ++blockIt)
   {
      const ScalarField_T* const fillField =
         blockIt->template getData< const ScalarField_T >(geometryHandler.fillFieldID_);

      // check if two ghost layers are required for the fill level field
      if (isObstacleInGlobalGhostLayer && fillField->nrOfGhostLayers() < uint_c(2) && geometryHandler.enableWetting_)
      {
         WALBERLA_ABORT(
            "With wetting enabled, the curvature computation with the finite difference method requires two ghost "
            "layers in the fill level field whenever solid obstacles are located in a global outermost ghost layer. "
            "For more information, see the remark in the description of ObstacleFillLevelSweep.h");
      }

      // get layout of fill level field (to be used in allocating the smoothed fill level field; cloning would
      // waste memory, as the fill level field might have two ghost layers, whereas the smoothed fill level field needs
      // only one ghost layer)
      fillFieldLayout = fillField->layout();
   }

   // IMPORTANT REMARK: ObstacleNormalSweep and ObstacleFillLevelSweep must be executed on all blocks, because the
   // SmoothingSweep requires meaningful values in the ghost layers.

   // add field for smoothed fill levels
   BlockDataID smoothFillFieldID = field::addToStorage< ScalarField_T >(
      geometryHandler.blockForest_, "Smooth fill level field", real_c(0), fillFieldLayout, uint_c(1));

   if (geometryHandler.enableWetting_)
   {
      // compute obstacle normals
      ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T > obstacleNormalSweep(
         geometryHandler.obstacleNormalFieldID_, geometryHandler.flagFieldID_, flagIDs::interfaceFlagID,
         flagIDs::liquidInterfaceGasFlagIDs, geometryHandler.obstacleFlagIDSet_, false, true, true);
      timeloop.add() << Sweep(obstacleNormalSweep, "Sweep: obstacle normal computation")
                     << AfterFunction(
                           Communication_T(geometryHandler.blockForest_, geometryHandler.obstacleNormalFieldID_),
                           "Communication: after obstacle normal sweep");

      // reflect fill level into obstacle cells such that they can be used for smoothing the fill level field and for
      // computing the interface normal; MUST be performed BEFORE SmoothingSweep
      ObstacleFillLevelSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > obstacleFillLevelSweep(
         smoothFillFieldID, geometryHandler.fillFieldID_, geometryHandler.flagFieldID_,
         geometryHandler.obstacleNormalFieldID_, flagIDs::liquidInterfaceGasFlagIDs,
         geometryHandler.obstacleFlagIDSet_);
      timeloop.add() << Sweep(obstacleFillLevelSweep, "Sweep: obstacle fill level computation")
                     << AfterFunction(Communication_T(geometryHandler.blockForest_, smoothFillFieldID),
                                      "Communication: after obstacle fill level sweep");
   }

   // smooth fill level field for decreasing error in finite difference normal and curvature computation (see
   // dissertation of S. Bogner, 2017 (section 4.4.2.1))
   SmoothingSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > smoothingSweep(
      smoothFillFieldID, geometryHandler.fillFieldID_, geometryHandler.flagFieldID_, flagIDs::liquidInterfaceGasFlagIDs,
      geometryHandler.obstacleFlagIDSet_, geometryHandler.enableWetting_);
   // IMPORTANT REMARK: SmoothingSweep must be executed on all blocks, because the algorithm works on all liquid,
   // interface and gas cells. This is necessary since the normals are not only computed in interface cells, but also in
   // the neighborhood of interface cells. Therefore, meaningful values for the fill levels of the second neighbors of
   // interface cells are also required in NormalSweep.
   timeloop.add() << Sweep(smoothingSweep, "Sweep: fill level smoothing")
                  << AfterFunction(Communication_T(geometryHandler.blockForest_, smoothFillFieldID),
                                   "Communication: after smoothing sweep");

   // compute interface normals (using smoothed fill level field)
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalSweep(
      geometryHandler.normalFieldID_, smoothFillFieldID, geometryHandler.flagFieldID_, flagIDs::interfaceFlagID,
      flagIDs::liquidInterfaceGasFlagIDs, geometryHandler.obstacleFlagIDSet_, true, geometryHandler.enableWetting_,
      true, geometryHandler.enableWetting_);
   timeloop.add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: normal")
                  << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.normalFieldID_),
                                   "Communication: after normal sweep");

   if (geometryHandler.computeCurvature_)
   {
      // compute interface curvature using finite differences according to Brackbill et al.
      CurvatureSweepFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvSweep(
         geometryHandler.curvatureFieldID_, geometryHandler.normalFieldID_, geometryHandler.obstacleNormalFieldID_,
         geometryHandler.flagFieldID_, flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs,
         geometryHandler.obstacleFlagIDSet_, geometryHandler.enableWetting_, geometryHandler.contactAngle_);
      timeloop.add() << Sweep(curvSweep, "Sweep: curvature computation (finite difference method)",
                              StateSweep::fullFreeSurface)
                     << Sweep(emptySweep(), "Empty sweep: curvature")
                     << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.curvatureFieldID_),
                                      "Communication: after curvature sweep");
   }
}

template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
void LocalTriangulation< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::addSweeps(
   SweepTimeloop& timeloop, const LocalTriangulation::SurfaceGeometryHandler_T& geometryHandler)
{
   using Communication_T = typename SurfaceGeometryHandler_T::Communication_T;
   using StateSweep      = typename SurfaceGeometryHandler_T::StateSweep;

   // compute interface normals
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalSweep(
      geometryHandler.normalFieldID_, geometryHandler.fillFieldID_, geometryHandler.flagFieldID_,
      flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs, geometryHandler.obstacleFlagIDSet_, false, false,
      true, false);
   timeloop.add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: normal")
                  << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.normalFieldID_),
                                   "Communication: after normal sweep");

   // compute obstacle normals
   ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T > obstacleNormalSweep(
      geometryHandler.obstacleNormalFieldID_, geometryHandler.flagFieldID_, flagIDs::interfaceFlagID,
      flagIDs::liquidInterfaceGasFlagIDs, geometryHandler.obstacleFlagIDSet_, true, false, false);
   timeloop.add() << Sweep(obstacleNormalSweep, "Sweep: obstacle normal computation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: obstacle normal")
                  << AfterFunction(
                        Communication_T(geometryHandler.blockForest_, geometryHandler.obstacleNormalFieldID_),
                        "Communication: after obstacle normal sweep");

   if (geometryHandler.computeCurvature_)
   {
      // compute interface curvature using local triangulation according to dissertation of T. Pohl, 2008
      CurvatureSweepLocalTriangulation< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvSweep(
         geometryHandler.blockForest_, geometryHandler.curvatureFieldID_, geometryHandler.normalFieldID_,
         geometryHandler.fillFieldID_, geometryHandler.flagFieldID_, geometryHandler.obstacleNormalFieldID_,
         flagIDs::interfaceFlagID, geometryHandler.obstacleFlagIDSet_, geometryHandler.enableWetting_,
         geometryHandler.contactAngle_);
      timeloop.add() << Sweep(curvSweep, "Sweep: curvature computation (local triangulation)",
                              StateSweep::fullFreeSurface)
                     << Sweep(emptySweep(), "Empty sweep: curvature")
                     << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.curvatureFieldID_),
                                      "Communication: after curvature sweep");
   }

   // sweep for detecting cells that need to be converted to interface cells for continuing the wetting
   // surface correctly
   // IMPORTANT REMARK: this MUST NOT be performed when using finite differences for curvature computation and can
   // otherwise lead to instabilities and errors
   if (geometryHandler.enableWetting_)
   {
      const auto& flagInfo = geometryHandler.freeSurfaceBoundaryHandling_->getFlagInfo();

      DetectWettingSweep<
         Stencil_T,
         typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::BoundaryHandling_T,
         FlagField_T, ScalarField_T, VectorField_T >
         detWetSweep(geometryHandler.freeSurfaceBoundaryHandling_->getHandlingID(), flagInfo,
                     geometryHandler.normalFieldID_, geometryHandler.fillFieldID_);
      timeloop.add() << Sweep(detWetSweep, "Sweep: wetting detection", StateSweep::fullFreeSurface)
                     << Sweep(emptySweep(), "Empty sweep: wetting detection")
                     << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.flagFieldID_),
                                      "Communication: after wetting detection sweep")
                     << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(geometryHandler.blockForest_,
                                                                                         geometryHandler.flagFieldID_),
                                      "Second ghost layer update: after wetting detection sweep (flag field)");
   }
}

template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
void SimpleFiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::addSweeps(
   SweepTimeloop& timeloop, const SimpleFiniteDifferenceMethod::SurfaceGeometryHandler_T& geometryHandler)
{
   using Communication_T = typename SurfaceGeometryHandler_T::Communication_T;
   using StateSweep      = typename SurfaceGeometryHandler_T::StateSweep;

   // compute interface normals
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalSweep(
      geometryHandler.normalFieldID_, geometryHandler.fillFieldID_, geometryHandler.flagFieldID_,
      flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs, geometryHandler.obstacleFlagIDSet_, false, false,
      false, false);
   timeloop.add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: normal")
                  << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.normalFieldID_),
                                   "Communication: after normal sweep");

   // extrapolation of normals to interface neighboring cells (required for computing the curvature with finite
   // differences)
   ExtrapolateNormalsSweep< Stencil_T, FlagField_T, VectorField_T > extNormalsSweep(
      geometryHandler.normalFieldID_, geometryHandler.flagFieldID_, flagIDs::interfaceFlagID);
   timeloop.add() << Sweep(extNormalsSweep, "Sweep: normal extrapolation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: normal extrapolation")
                  << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.normalFieldID_),
                                   "Communication: after normal extrapolation sweep");

   if (geometryHandler.computeCurvature_)
   {
      // curvature computation using finite differences
      CurvatureSweepSimpleFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvSweep(
         geometryHandler.curvatureFieldID_, geometryHandler.normalFieldID_, geometryHandler.flagFieldID_,
         flagIDs::interfaceFlagID, geometryHandler.obstacleFlagIDSet_, geometryHandler.enableWetting_,
         geometryHandler.contactAngle_);
      timeloop.add() << Sweep(curvSweep, "Sweep: curvature computation (simple finite difference method)",
                              StateSweep::fullFreeSurface)
                     << Sweep(emptySweep(), "Empty sweep: curvature")
                     << AfterFunction(Communication_T(geometryHandler.blockForest_, geometryHandler.curvatureFieldID_),
                                      "Communication: after curvature sweep ");
   }
}

} // namespace curvature_model
} // namespace free_surface
} // namespace walberla