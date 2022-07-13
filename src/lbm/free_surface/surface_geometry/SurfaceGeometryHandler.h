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
//! \file SurfaceGeometryHandler.h
//! \ingroup surface_geometry
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Handles the surface geometry (normal and curvature computation) by creating fields and adding sweeps.
//
//======================================================================================================================

#pragma once

#include "core/StringUtility.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/FlagInfo.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <type_traits>
#include <vector>

#include "CurvatureModel.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Handles the surface geometry (normal and curvature computation) by creating fields and adding sweeps.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class SurfaceGeometryHandler
{
 protected:
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;
   using vector_t                      = typename std::remove_const< typename VectorField_T::value_type >::type;

   // explicitly use either D2Q9 or D3Q27 here, as the geometry operations require (or are most accurate with) the full
   // neighborhood;
   using Stencil_T =
      typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
   using Communication_T = blockforest::SimpleCommunication< Stencil_T >;
   using StateSweep      = BlockStateDetectorSweep< FlagField_T >; // used in friend classes

 public:
   SurfaceGeometryHandler(const std::shared_ptr< StructuredBlockForest >& blockForest,
                          const std::shared_ptr< FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                          const BlockDataID& fillFieldID, const std::string& curvatureModel, bool computeCurvature,
                          bool enableWetting, real_t contactAngleInDegrees)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), fillFieldID_(fillFieldID),
        curvatureModel_(curvatureModel), computeCurvature_(computeCurvature), enableWetting_(enableWetting),
        contactAngle_(ContactAngle(contactAngleInDegrees))
   {
      curvatureFieldID_ =
         field::addToStorage< ScalarField_T >(blockForest_, "Curvature field", real_c(0), field::fzyx, uint_c(1));
      normalFieldID_         = field::addToStorage< VectorField_T >(blockForest_, "Normal field", vector_t(real_c(0)),
                                                            field::fzyx, uint_c(1));
      obstacleNormalFieldID_ = field::addToStorage< VectorField_T >(blockForest_, "Obstacle normal field",
                                                                    vector_t(real_c(0)), field::fzyx, uint_c(1));

      flagFieldID_ = freeSurfaceBoundaryHandling_->getFlagFieldID();

      obstacleFlagIDSet_ = freeSurfaceBoundaryHandling_->getFlagInfo().getObstacleIDSet();

      if (LatticeModel_T::Stencil::D == uint_t(2))
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "IMPORTANT REMARK: You are using a D2Q9 stencil in SurfaceGeometryHandler. Be aware that the "
            "results might slightly differ when compared to a D3Q19 stencil and periodicity in the third direction. "
            "This is caused by the smoothing of the fill level field, where the additional directions in the D3Q27 add "
            "additional weights to the smoothing kernel. Therefore, the resulting smoothed fill level will be "
            "different.")
      }
   }

   ConstBlockDataID getConstCurvatureFieldID() const { return curvatureFieldID_; }
   ConstBlockDataID getConstNormalFieldID() const { return normalFieldID_; }
   ConstBlockDataID getConstObstNormalFieldID() const { return obstacleNormalFieldID_; }

   BlockDataID getCurvatureFieldID() const { return curvatureFieldID_; }
   BlockDataID getNormalFieldID() const { return normalFieldID_; }
   BlockDataID getObstNormalFieldID() const { return obstacleNormalFieldID_; }

   void addSweeps(SweepTimeloop& timeloop) const
   {
      auto blockStateUpdate = StateSweep(blockForest_, freeSurfaceBoundaryHandling_->getFlagInfo(), flagFieldID_);

      if (!string_icompare(curvatureModel_, "FiniteDifferenceMethod"))
      {
         curvature_model::FiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            model;
         model.addSweeps(timeloop, *this);
      }
      else
      {
         if (!string_icompare(curvatureModel_, "LocalTriangulation"))
         {
            curvature_model::LocalTriangulation< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
               model;
            model.addSweeps(timeloop, *this);
         }
         else
         {
            if (!string_icompare(curvatureModel_, "SimpleFiniteDifferenceMethod"))
            {
               curvature_model::SimpleFiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T,
                                                              VectorField_T >
                  model;
               model.addSweeps(timeloop, *this);
            }
            else { WALBERLA_ABORT("The specified curvature model is unknown.") }
         }
      }
   }

 protected:
   std::shared_ptr< StructuredBlockForest > blockForest_; // used by friend classes

   std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   BlockDataID flagFieldID_;
   ConstBlockDataID fillFieldID_;

   BlockDataID curvatureFieldID_;
   BlockDataID normalFieldID_;
   BlockDataID obstacleNormalFieldID_; // mean normal in obstacle cells required for e.g. artificial curvature contact
                                       // model (dissertation of S. Donath, 2011, section 6.5.3.2)

   Set< FlagUID > obstacleFlagIDSet_; // used by friend classes (see CurvatureModel.impl.h)

   std::string curvatureModel_;
   bool computeCurvature_;     // allows to not compute curvature (just normal) when e.g. the surface tension is 0
   bool enableWetting_;        // used by friend classes
   ContactAngle contactAngle_; // used by friend classes

   friend class curvature_model::FiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T,
                                                         VectorField_T >;
   friend class curvature_model::LocalTriangulation< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T,
                                                     VectorField_T >;
   friend class curvature_model::SimpleFiniteDifferenceMethod< Stencil_T, LatticeModel_T, FlagField_T, ScalarField_T,
                                                               VectorField_T >;
}; // class SurfaceGeometryHandler

} // namespace free_surface
} // namespace walberla
