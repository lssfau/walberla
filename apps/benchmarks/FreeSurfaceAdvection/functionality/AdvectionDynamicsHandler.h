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
//! \file AdvectionDynamicsHandler.h
//! \ingroup surface_dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Handles free surface advection (without LBM flow simulation; this is a simplified SurfaceDynamicsHandler).
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
#include "lbm/free_surface/dynamics/CellConversionSweep.h"
#include "lbm/free_surface/dynamics/ConversionFlagsResetSweep.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionModel.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionSweep.h"
#include "lbm/free_surface/dynamics/PdfReconstructionModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "AdvectSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class AdvectionDynamicsHandler
{
 protected:
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::Stencil >;

   // communication in corner directions (D2Q9/D3Q27) is required for all fields but the PDF field
   using CommunicationStencil_T =
      typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
   using CommunicationCorner_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

 public:
   AdvectionDynamicsHandler(const std::shared_ptr< StructuredBlockForest >& blockForest, BlockDataID pdfFieldID,
                            BlockDataID flagFieldID, BlockDataID fillFieldID, ConstBlockDataID normalFieldID,
                            const std::shared_ptr< FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                            const std::shared_ptr< BubbleModelBase >& bubbleModel,
                            const std::string& pdfReconstructionModel, const std::string& excessMassDistributionModel,
                            bool useSimpleMassExchange, real_t cellConversionThreshold,
                            real_t cellConversionForceThreshold)
      : blockForest_(blockForest), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID),
        normalFieldID_(normalFieldID), bubbleModel_(bubbleModel),
        freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfReconstructionModel_(pdfReconstructionModel),
        excessMassDistributionModel_({ excessMassDistributionModel }), useSimpleMassExchange_(useSimpleMassExchange),
        cellConversionThreshold_(cellConversionThreshold), cellConversionForceThreshold_(cellConversionForceThreshold)
   {
      WALBERLA_CHECK(LatticeModel_T::compressible,
                     "The free surface lattice Boltzmann extension works only with compressible LBM models.");

      if (excessMassDistributionModel_.isEvenlyAllInterfaceFallbackLiquidType())
      {
         // add additional field for storing excess mass in liquid cells
         excessMassFieldID_ =
            field::addToStorage< ScalarField_T >(blockForest_, "Excess mass", real_c(0), field::fzyx, uint_c(1));
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

      // sweep for
      // - reconstruction of PDFs in interface cells
      // - streaming of PDFs in interface cells (and liquid cells on the same block)
      // - advection of mass
      // - update bubble volumes
      // - marking interface cells for conversion
      const AdvectSweep< LatticeModel_T, typename FreeSurfaceBoundaryHandling_T::BoundaryHandling_T, FlagField_T,
                         typename FreeSurfaceBoundaryHandling_T::FlagInfo_T, ScalarField_T, VectorField_T >
         advectSweep(freeSurfaceBoundaryHandling_->getHandlingID(), fillFieldID_, flagFieldID_, pdfFieldID_, flagInfo,
                     pdfReconstructionModel_, useSimpleMassExchange_, cellConversionThreshold_,
                     cellConversionForceThreshold_);
      // sweep acts only on blocks with at least one interface cell (due to StateSweep::fullFreeSurface)
      timeloop.add()
         << Sweep(advectSweep, "Sweep: Advect", StateSweep::fullFreeSurface)
         << Sweep(emptySweep, "Empty sweep: Advect")
         // do not communicate PDFs here:
         // - stream on blocks with "StateSweep::fullFreeSurface" was performed here using post-collision PDFs
         // - stream on other blocks is performed below and should also use post-collision PDFs
         // => if PDFs were communicated here, the ghost layer of other blocks would have post-stream PDFs
         << AfterFunction(CommunicationCorner_T(blockForest_, fillFieldID_, flagFieldID_),
                          "Communication: after Advect sweep")
         << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest_, fillFieldID_),
                          "Second ghost layer update: after Advect sweep (fill level field)")
         << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest_, flagFieldID_),
                          "Second ghost layer update: after Advect sweep (flag field)");

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

      // distribute excess mass:
      // - excess mass: mass that is free after conversion from interface to gas/liquid cells
      // - update the bubble model
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
            if (excessMassDistributionModel_.isEvenlyAllInterfaceFallbackLiquidType())
            {
               const ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T,
                                                                    VectorField_T >
                  distributeMassSweep(excessMassDistributionModel_, fillFieldID_, flagFieldID_, pdfFieldID_, flagInfo,
                                      excessMassFieldID_);
               timeloop.add()
                  // perform this sweep also on "onlyLBM" blocks because liquid cells also exchange excess mass here
                  << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                  << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::onlyLBM)
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

   ConstBlockDataID normalFieldID_;

   std::shared_ptr< BubbleModelBase > bubbleModel_;
   std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   PdfReconstructionModel pdfReconstructionModel_;
   ExcessMassDistributionModel excessMassDistributionModel_;

   bool useSimpleMassExchange_;
   real_t cellConversionThreshold_;
   real_t cellConversionForceThreshold_;

   BlockDataID excessMassFieldID_ = BlockDataID();
}; // class AdvectionDynamicsHandler

} // namespace free_surface
} // namespace walberla