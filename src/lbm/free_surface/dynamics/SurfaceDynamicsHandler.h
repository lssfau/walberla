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
//! \file SurfaceDynamicsHandler.h
//! \ingroup surface_dynamics
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Handles the free surface dynamics (mass advection, LBM, boundary condition, cell conversion etc.).
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/blockforest/communication/UpdateSecondGhostLayer.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"
#include "lbm/lattice_model/SmagorinskyLES.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "CellConversionSweep.h"
#include "ConversionFlagsResetSweep.h"
#include "ExcessMassDistributionModel.h"
#include "ExcessMassDistributionSweep.h"
#include "ForceWeightingSweep.h"
#include "PdfReconstructionModel.h"
#include "PdfRefillingModel.h"
#include "PdfRefillingSweep.h"
#include "StreamReconstructAdvectSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T,
          bool useCodegen = false >
class SurfaceDynamicsHandler
{
 protected:
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::Stencil >;

   // communication in corner directions (D2Q9/D3Q27) is required for all fields but the PDF field
   using CommunicationStencil_T =
      typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
   using CommunicationCorner_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

 public:
   SurfaceDynamicsHandler(const std::shared_ptr< StructuredBlockForest >& blockForest, BlockDataID pdfFieldID,
                          BlockDataID flagFieldID, BlockDataID fillFieldID, BlockDataID forceFieldID,
                          ConstBlockDataID normalFieldID, ConstBlockDataID curvatureFieldID,
                          const std::shared_ptr< FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                          const std::shared_ptr< BubbleModelBase >& bubbleModel,
                          const std::string& pdfReconstructionModel, const std::string& pdfRefillingModel,
                          const std::string& excessMassDistributionModel, real_t relaxationRate,
                          const Vector3< real_t >& globalForce, real_t surfaceTension, bool enableForceWeighting,
                          bool useSimpleMassExchange, real_t cellConversionThreshold,
                          real_t cellConversionForceThreshold, BlockDataID relaxationRateFieldID = BlockDataID(),
                          real_t smagorinskyConstant = real_c(0))
      : blockForest_(blockForest), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID),
        forceFieldID_(forceFieldID), normalFieldID_(normalFieldID), curvatureFieldID_(curvatureFieldID),
        bubbleModel_(bubbleModel), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling),
        pdfReconstructionModel_(pdfReconstructionModel), pdfRefillingModel_({ pdfRefillingModel }),
        excessMassDistributionModel_({ excessMassDistributionModel }), relaxationRate_(relaxationRate),
        globalForce_(globalForce), surfaceTension_(surfaceTension), enableForceWeighting_(enableForceWeighting),
        useSimpleMassExchange_(useSimpleMassExchange), cellConversionThreshold_(cellConversionThreshold),
        cellConversionForceThreshold_(cellConversionForceThreshold), relaxationRateFieldID_(relaxationRateFieldID),
        smagorinskyConstant_(smagorinskyConstant)
   {
      WALBERLA_CHECK(LatticeModel_T::compressible,
                     "The free surface lattice Boltzmann extension works only with compressible LBM models.");

      if (useCodegen && !realIsEqual(smagorinskyConstant_, real_c(0)))
      {
         WALBERLA_ABORT("When using a generated LBM kernel from lbmpy, please use lbmpy's inbuilt-functionality to "
                        "generate the Smagorinsky model directly into the kernel.");
      }

      if (excessMassDistributionModel_.isEvenlyLiquidAndAllInterfacePreferInterfaceType())
      {
         // add additional field for storing excess mass in liquid cells
         excessMassFieldID_ =
            field::addToStorage< ScalarField_T >(blockForest_, "Excess mass", real_c(0), field::fzyx, uint_c(1));
      }

      if (LatticeModel_T::Stencil::D == uint_t(2))
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "IMPORTANT REMARK: You are using a D2Q9 stencil in SurfaceDynamicsHandler. In the FSLBM, a D2Q9 setup is "
            "not identical to a D3Q19 setup with periodicity in the third direction but only identical to the same "
            "D3Q27 setup. This comes from the distribution of excess mass, where the number of available neighbors "
            "differs for D2Q9 and a periodic D3Q19 setup.")
      }
   }

   ConstBlockDataID getConstExcessMassFieldID() const { return excessMassFieldID_; }

   /********************************************************************************************************************
    * The order of these sweeps is similar to page 40 in the dissertation of T. Pohl, 2008.
    *******************************************************************************************************************/
   void addSweeps(SweepTimeloop& timeloop) const
   {
      using StateSweep = BlockStateDetectorSweep< FlagField_T >;

      const auto& flagInfo = freeSurfaceBoundaryHandling_->getFlagInfo();

      const auto blockStateUpdate = StateSweep(blockForest_, flagInfo, flagFieldID_);

      // empty sweeps required for using selectors (e.g. StateSweep::onlyGasAndBoundary)
      const auto emptySweep = [](IBlock*) {};

      // add standard waLBerla boundary handling
      timeloop.add() << Sweep(freeSurfaceBoundaryHandling_->getBoundarySweep(), "Sweep: boundary handling",
                              Set< SUID >::emptySet(), StateSweep::onlyGasAndBoundary)
                     << Sweep(emptySweep, "Empty sweep: boundary handling", StateSweep::onlyGasAndBoundary);

      if (enableForceWeighting_)
      {
         // add sweep for weighting force in interface cells with fill level and density
         const ForceWeightingSweep< LatticeModel_T, FlagField_T, VectorField_T, ScalarField_T > forceWeightingSweep(
            forceFieldID_, pdfFieldID_, flagFieldID_, fillFieldID_, flagInfo, globalForce_);
         timeloop.add() << Sweep(forceWeightingSweep, "Sweep: force weighting", Set< SUID >::emptySet(),
                                 StateSweep::onlyGasAndBoundary)
                        << Sweep(emptySweep, "Empty sweep: force weighting", StateSweep::onlyGasAndBoundary)
                        << AfterFunction(CommunicationCorner_T(blockForest_, forceFieldID_),
                                         "Communication: after force weighting sweep");
      }

      // sweep for
      // - reconstruction of PDFs in interface cells
      // - streaming of PDFs in interface cells (and liquid cells on the same block)
      // - advection of mass
      // - update bubble volumes
      // - marking interface cells for conversion
      const StreamReconstructAdvectSweep< LatticeModel_T, typename FreeSurfaceBoundaryHandling_T::BoundaryHandling_T,
                                          FlagField_T, typename FreeSurfaceBoundaryHandling_T::FlagInfo_T,
                                          ScalarField_T, VectorField_T, useCodegen >
         streamReconstructAdvectSweep(surfaceTension_, freeSurfaceBoundaryHandling_->getHandlingID(), fillFieldID_,
                                      flagFieldID_, pdfFieldID_, normalFieldID_, curvatureFieldID_, flagInfo,
                                      bubbleModel_.get(), pdfReconstructionModel_, useSimpleMassExchange_,
                                      cellConversionThreshold_, cellConversionForceThreshold_);
      // sweep acts only on blocks with at least one interface cell (due to StateSweep::fullFreeSurface)
      timeloop.add()
         << Sweep(streamReconstructAdvectSweep, "Sweep: StreamReconstructAdvect", StateSweep::fullFreeSurface)
         << Sweep(emptySweep, "Empty sweep: StreamReconstructAdvect")
         // do not communicate PDFs here:
         // - stream on blocks with "StateSweep::fullFreeSurface" was performed here using post-collision PDFs
         // - stream on other blocks is performed below and should also use post-collision PDFs
         // => if PDFs were communicated here, the ghost layer of other blocks would have post-stream PDFs
         << AfterFunction(CommunicationCorner_T(blockForest_, fillFieldID_, flagFieldID_),
                          "Communication: after StreamReconstructAdvect sweep")
         << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest_, fillFieldID_),
                          "Second ghost layer update: after StreamReconstructAdvect sweep (fill level field)")
         << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest_, flagFieldID_),
                          "Second ghost layer update: after StreamReconstructAdvect sweep (flag field)");

      if constexpr (useCodegen)
      {
         auto lbmSweepGenerated = typename LatticeModel_T::Sweep(pdfFieldID_);

         // temporary class for being able to call the LBM collision with operator()
         class CollideSweep
         {
          public:
            CollideSweep(const typename LatticeModel_T::Sweep& sweep) : sweep_(sweep){};

            void operator()(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0))
            {
               sweep_.collide(block, numberOfGhostLayersToInclude);
            }

          private:
            typename LatticeModel_T::Sweep sweep_;
         };

         timeloop.add() << Sweep(CollideSweep(lbmSweepGenerated), "Sweep: collision (generated)",
                                 StateSweep::fullFreeSurface)
                        << Sweep(lbmSweepGenerated, "Sweep: streamCollide (generated)", StateSweep::onlyLBM)
                        << Sweep(emptySweep, "Empty sweep: streamCollide (generated)")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after streamCollide (generated)");
      }
      else
      {
         // sweep for standard LBM stream and collision
         const auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >(pdfFieldID_, flagFieldID_,
                                                                                     flagIDs::liquidInterfaceFlagIDs);

         // LBM collision in interface cells; standard LBM stream and collision in liquid cells
         if (!realIsEqual(smagorinskyConstant_, real_c(0), real_c(1e-14))) // using Smagorinsky turbulence model
         {
            const real_t kinematicViscosity = (real_c(1) / relaxationRate_ - real_c(0.5)) / real_c(3);

            // standard LBM stream in liquid cells that have not been streamed, yet
            timeloop.add() << Sweep(lbm::makeStreamSweep(lbmSweep), "Stream sweep", StateSweep::onlyLBM)
                           << Sweep(emptySweep, "Deactivated Stream sweep")
                           << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                            "Communication after Stream sweep");

            // sweep for turbulence modelling
            const lbm::SmagorinskyLES< LatticeModel_T > smagorinskySweep(
               blockForest_, pdfFieldID_, relaxationRateFieldID_, kinematicViscosity, real_c(0.12));

            timeloop.add()
               // Smagorinsky turbulence model
               << BeforeFunction(smagorinskySweep, "Sweep: Smagorinsky turbulence model")
               << BeforeFunction(CommunicationCorner_T(blockForest_, relaxationRateFieldID_),
                                 "Communication: after Smagorinsky sweep")
               // standard LBM collision
               << Sweep(lbm::makeCollideSweep(lbmSweep), "Sweep: collision after Smagorinsky sweep")
               << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                "Communication: after collision sweep with preceding Smagorinsky sweep");
         }
         else
         {
            // no turbulence model
            timeloop.add()
               // collision in interface cells and liquid cells that have already been streamed (in
               // streamReconstructAdvectSweep due to StateSweep::fullFreeSurface)
               << Sweep(lbm::makeCollideSweep(lbmSweep), "Sweep: collision", StateSweep::fullFreeSurface)
               // standard LBM stream-collide in liquid cells that have not been streamed, yet
               << Sweep(*lbmSweep, "Sweep: streamCollide", StateSweep::onlyLBM)
               << Sweep(emptySweep, "Empty sweep: streamCollide")
               << AfterFunction(Communication_T(blockForest_, pdfFieldID_), "Communication: after streamCollide sweep");
         }
      }

      // convert cells
      // - according to the flags from StreamReconstructAdvectSweep (interface -> gas/liquid)
      // - to ensure a closed layer of interface cells (gas/liquid -> interface)
      // - detect and register bubble merges/splits (bubble volumes are already updated in StreamReconstructAdvectSweep)
      // - convert cells and initialize PDFs near inflow boundaries
      const CellConversionSweep< LatticeModel_T, typename FreeSurfaceBoundaryHandling_T::BoundaryHandling_T,
                                 ScalarField_T >
         cellConvSweep(freeSurfaceBoundaryHandling_->getHandlingID(), pdfFieldID_, flagInfo, bubbleModel_.get());
      timeloop.add() << Sweep(cellConvSweep, "Sweep: cell conversion", StateSweep::fullFreeSurface)
                     << Sweep(emptySweep, "Empty sweep: cell conversion")
                     << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                      "Communication: after cell conversion sweep (PDF field)")
                     // communicate the flag field also in corner directions
                     << AfterFunction(CommunicationCorner_T(blockForest_, flagFieldID_),
                                      "Communication: after cell conversion sweep (flag field)")
                     << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest_, flagFieldID_),
                                      "Second ghost layer update: after cell conversion sweep (flag field)");

      // reinitialize PDFs, i.e., refill cells that were converted from gas to interface
      // - when the flag "convertedFromGasToInterface" has been set (by CellConversionSweep)
      // - according to the method specified with pdfRefillingModel_
      switch (pdfRefillingModel_.getModelType())
      { // the scope for each "case" is required since variables are defined within "case"
      case PdfRefillingModel::RefillingModel::EquilibriumRefilling: {
         const EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > equilibriumRefillingSweep(
            pdfFieldID_, flagFieldID_, flagInfo, true);
         timeloop.add() << Sweep(equilibriumRefillingSweep, "Sweep: EquilibriumRefilling", StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: EquilibriumRefilling")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after EquilibriumRefilling sweep");
         break;
      }
      case PdfRefillingModel::RefillingModel::AverageRefilling: {
         const AverageRefillingSweep< LatticeModel_T, FlagField_T > averageRefillingSweep(pdfFieldID_, flagFieldID_,
                                                                                          flagInfo, true);
         timeloop.add() << Sweep(averageRefillingSweep, "Sweep: AverageRefilling", StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: AverageRefilling")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after AverageRefilling sweep");
         break;
      }
      case PdfRefillingModel::RefillingModel::EquilibriumAndNonEquilibriumRefilling: {
         // default: extrapolation order: 0
         const EquilibriumAndNonEquilibriumRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            equilibriumAndNonEquilibriumRefillingSweep(pdfFieldID_, flagFieldID_, fillFieldID_, flagInfo, uint_c(0),
                                                       true);
         timeloop.add() << Sweep(equilibriumAndNonEquilibriumRefillingSweep,
                                 "Sweep: EquilibriumAndNonEquilibriumRefilling sweep", StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: EquilibriumAndNonEquilibriumRefilling")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after EquilibriumAndNonEquilibriumRefilling sweep");
         break;
      }
      case PdfRefillingModel::RefillingModel::ExtrapolationRefilling: {
         // default: extrapolation order: 2
         const ExtrapolationRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            extrapolationRefillingSweep(pdfFieldID_, flagFieldID_, fillFieldID_, flagInfo, uint_c(2), true);
         timeloop.add() << Sweep(extrapolationRefillingSweep, "Sweep: ExtrapolationRefilling",
                                 StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: ExtrapolationRefilling")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after ExtrapolationRefilling sweep");
         break;
      }
      case PdfRefillingModel::RefillingModel::GradsMomentsRefilling: {
         const GradsMomentsRefillingSweep< LatticeModel_T, FlagField_T > gradsMomentRefillingSweep(
            pdfFieldID_, flagFieldID_, flagInfo, relaxationRate_, true);
         timeloop.add() << Sweep(gradsMomentRefillingSweep, "Sweep: GradsMomentRefilling", StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: GradsMomentRefilling")
                        << AfterFunction(Communication_T(blockForest_, pdfFieldID_),
                                         "Communication: after GradsMomentRefilling sweep");
         break;
      }
      default:
         WALBERLA_ABORT("The specified pdf refilling model is not available.");
      }

      // distribute excess mass:
      // - excess mass: mass that is free after conversion from interface to gas/liquid cells
      // - update the bubble model
      // IMPORTANT REMARK: this sweep computes the mass via the density, i.e., the PDF field must be up-to-date and the
      // PdfRefillingSweep must have been performed
      if (excessMassDistributionModel_.isEvenlyType())
      {
         const ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            distributeMassSweep(excessMassDistributionModel_, fillFieldID_, flagFieldID_, pdfFieldID_, flagInfo);
         timeloop.add() << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                        << Sweep(emptySweep, "Empty sweep: distribute excess mass")
                        << AfterFunction(CommunicationCorner_T(blockForest_, fillFieldID_),
                                         "Communication: after excess mass distribution sweep")
                        << AfterFunction(
                              blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest_, fillFieldID_),
                              "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                        // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits are
                        // already detected and registered by CellConversionSweep
                        << AfterFunction(std::bind(&bubble_model::BubbleModelBase::update, bubbleModel_),
                                         "Sweep: bubble model update");
      }
      else
      {
         if (excessMassDistributionModel_.isWeightedType())
         {
            const ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T,
                                                                VectorField_T >
               distributeMassSweep(excessMassDistributionModel_, fillFieldID_, flagFieldID_, pdfFieldID_, flagInfo,
                                   normalFieldID_);
            timeloop.add() << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                           << Sweep(emptySweep, "Empty sweep: distribute excess mass")
                           << AfterFunction(CommunicationCorner_T(blockForest_, fillFieldID_),
                                            "Communication: after excess mass distribution sweep")
                           << AfterFunction(
                                 blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest_, fillFieldID_),
                                 "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                           // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits
                           // are already detected and registered by CellConversionSweep
                           << AfterFunction(std::bind(&bubble_model::BubbleModelBase::update, bubbleModel_),
                                            "Sweep: bubble model update");
         }
         else
         {
            if (excessMassDistributionModel_.isEvenlyLiquidAndAllInterfacePreferInterfaceType())
            {
               const ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T,
                                                                    VectorField_T >
                  distributeMassSweep(excessMassDistributionModel_, fillFieldID_, flagFieldID_, pdfFieldID_, flagInfo,
                                      excessMassFieldID_);
               timeloop.add()
                  << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep, "Empty sweep: distribute excess mass")
                  << AfterFunction(CommunicationCorner_T(blockForest_, fillFieldID_, excessMassFieldID_),
                                   "Communication: after excess mass distribution sweep")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest_, fillFieldID_),
                                   "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                  // update bubble model, i.e., perform registered bubble merges/splits; bubble
                  // merges/splits are already detected and registered by CellConversionSweep
                  << AfterFunction(std::bind(&bubble_model::BubbleModelBase::update, bubbleModel_),
                                   "Sweep: bubble model update");
            }
         }
      }

      // reset all flags that signal cell conversions (except "keepInterfaceForWettingFlag")
      ConversionFlagsResetSweep< FlagField_T > resetConversionFlagsSweep(flagFieldID_, flagInfo);
      timeloop.add() << Sweep(resetConversionFlagsSweep, "Sweep: conversion flag reset", StateSweep::fullFreeSurface)
                     << Sweep(emptySweep, "Empty sweep: conversion flag reset")
                     << AfterFunction(CommunicationCorner_T(blockForest_, flagFieldID_),
                                      "Communication: after excess mass distribution sweep")
                     << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest_, flagFieldID_),
                                      "Second ghost layer update: after excess mass distribution sweep (flag field)");

      // update block states
      timeloop.add() << Sweep(blockStateUpdate, "Sweep: block state update");
   }

 private:
   std::shared_ptr< StructuredBlockForest > blockForest_;

   BlockDataID pdfFieldID_;
   BlockDataID flagFieldID_;
   BlockDataID fillFieldID_;
   BlockDataID forceFieldID_;

   ConstBlockDataID normalFieldID_;
   ConstBlockDataID curvatureFieldID_;

   std::shared_ptr< BubbleModelBase > bubbleModel_;
   std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   PdfReconstructionModel pdfReconstructionModel_;
   PdfRefillingModel pdfRefillingModel_;
   ExcessMassDistributionModel excessMassDistributionModel_;
   real_t relaxationRate_;
   Vector3< real_t > globalForce_;
   real_t surfaceTension_;
   bool enableForceWeighting_;
   bool useSimpleMassExchange_;
   real_t cellConversionThreshold_;
   real_t cellConversionForceThreshold_;

   BlockDataID relaxationRateFieldID_;
   real_t smagorinskyConstant_;

   BlockDataID excessMassFieldID_ = BlockDataID();
}; // class SurfaceDynamicsHandler

} // namespace free_surface
} // namespace walberla