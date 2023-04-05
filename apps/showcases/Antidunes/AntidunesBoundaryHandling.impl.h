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
//! \file AntidunesBoundaryHandling.impl.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Boundary handling for the free surface LBM module.
//
//======================================================================================================================

#pragma once

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"

#include "geometry/initializer/BoundaryFromCellInterval.h"
#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "geometry/initializer/BoundaryFromImage.h"
#include "geometry/structured/GrayScaleImage.h"

#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/InterfaceFromFillLevel.h"
#include "lbm/lattice_model/CollisionModel.h"

#include "AntidunesBoundaryHandling.h"

namespace walberla
{
namespace antidunes
{
namespace free_surface
{
namespace internal
{
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
class BoundaryBlockDataHandling
   : public domain_decomposition::BlockDataHandling< typename AntidunesBoundaryHandling<
        LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::BoundaryHandling_T >
{
 public:
   using BoundaryHandling_T =
      typename AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                          ParticleAccessor_T >::BoundaryHandling_T; // handling as defined in
                                                                                    // AntidunesBoundaryHandling.h

   BoundaryBlockDataHandling(
      const AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >* boundary)
      : boundary_(boundary)
   {}

   // initialize standard waLBerla boundary handling
   BoundaryHandling_T* initialize(IBlock* const block)
   {
      using B      = AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >;
      using flag_t = typename B::flag_t;

      // get fields
      FlagField_T* const flagField           = block->getData< FlagField_T >(boundary_->getFlagFieldID());
      typename B::PdfField_T* const pdfField = block->getData< typename B::PdfField_T >(boundary_->getPdfFieldID());
      lbm_mesapd_coupling::ParticleField_T* const particleField =
         block->getData< lbm_mesapd_coupling::ParticleField_T >(boundary_->getParticleFieldID());

      auto interfaceFlag = flag_t(flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID));
      auto liquidFlag    = flag_t(flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID));

      // domainMask is used to identify liquid and interface cells
      auto domainMask = flag_t(liquidFlag | interfaceFlag);
      WALBERLA_ASSERT(domainMask != 0);

      // initialize boundary conditions
      typename B::UBB_T ubb(B::ubbBoundaryID, B::ubbFlagID, pdfField, flagField);
      typename B::UBB_Inflow_T ubbInflow(B::ubbInflowBoundaryID, B::ubbInflowFlagID, pdfField, flagField);
      typename B::NoSlip_T noslip(B::noSlipBoundaryID, B::noSlipFlagID, pdfField);
      typename B::Pressure_T pressure(B::pressureBoundaryID, B::pressureFlagID, block, pdfField, flagField,
                                      interfaceFlag, real_c(1.0));
      typename B::Pressure_T pressureOutflow(B::pressureOutflowBoundaryID, B::pressureOutflowFlagID, block, pdfField,
                                             flagField, interfaceFlag, real_c(1.0));
      typename B::Outlet_T outlet(B::outletBoundaryID, B::outletFlagID, pdfField, flagField, domainMask);
      typename B::FreeSlip_T freeSlip(B::freeSlipBoundaryID, B::freeSlipFlagID, pdfField, flagField, domainMask);
      typename B::MovingObstacle_T curvedLinear(B::movingObstacleBoundaryID, B::movingObstacleFlagID, pdfField,
                                                flagField, particleField, boundary_->getParticleAccessor(), domainMask,
                                                *boundary_->getBlockForest(), *block,
                                                boundary_->getHydrostaticDensityFct());

      return new BoundaryHandling_T("Boundary Handling", flagField, domainMask, noslip, ubb, ubbInflow, pressure,
                                    pressureOutflow, outlet, freeSlip, curvedLinear);
   }

   void serialize(IBlock* const block, const BlockDataID& id, mpi::SendBuffer& buffer)
   {
      BoundaryHandling_T* const boundaryHandlingPtr = block->getData< BoundaryHandling_T >(id);
      CellInterval everyCell                        = boundaryHandlingPtr->getFlagField()->xyzSizeWithGhostLayer();
      boundaryHandlingPtr->pack(buffer, everyCell, true);
   }

   BoundaryHandling_T* deserialize(IBlock* const block) { return initialize(block); }

   void deserialize(IBlock* const block, const BlockDataID& id, mpi::RecvBuffer& buffer)
   {
      BoundaryHandling_T* const boundaryHandlingPtr = block->getData< BoundaryHandling_T >(id);
      CellInterval everyCell                        = boundaryHandlingPtr->getFlagField()->xyzSizeWithGhostLayer();
      boundaryHandlingPtr->unpack(buffer, everyCell, true);
   }

 private:
   const AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >* boundary_;
}; // class BoundaryBlockDataHandling

// helper function wrapper for adding the flag field to the block storage; since the input parameter for an
// initialization function in field::addFlagFieldToStorage() is a std::function<void(FlagField_T*,IBlock* const)>, we
// need a function wrapper that has both these input parameters; as FlagInfo< FlagField_T >::registerFlags() does not
// have both of these input parameters, a function wrapper with both input parameters is created and the second input
// parameter is simply ignored inside the function wrapper
template< typename FlagField_T >
void flagFieldInitFunction(FlagField_T* flagField, IBlock* const, const Set< field::FlagUID >& obstacleIDs,
                           const Set< field::FlagUID >& outflowIDs, const Set< field::FlagUID >& inflowIDs,
                           const Set< field::FlagUID >& freeSlipIDs)
{
   // register flags in the flag field
   walberla::free_surface::FlagInfo< FlagField_T >::registerFlags(flagField, obstacleIDs, outflowIDs, inflowIDs,
                                                                  freeSlipIDs);
}

} // namespace internal

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::AntidunesBoundaryHandling(
   const std::shared_ptr< StructuredBlockForest >& blockForest, BlockDataID pdfFieldID, BlockDataID fillLevelID,
   BlockDataID particleFieldID, const shared_ptr< ParticleAccessor_T >& ac,
   std::function< real_t(const Vector3< real_t >&) > hydrostaticDensityFct)
   : blockForest_(blockForest), pdfFieldID_(pdfFieldID), fillFieldID_(fillLevelID), particleFieldID_(particleFieldID),
     ac_(ac), hydrostaticDensityFct_(std::move(hydrostaticDensityFct)), comm_(blockForest)
{
   // initialize obstacleIDs
   Set< FlagUID > obstacleIDs;
   obstacleIDs += noSlipFlagID;
   obstacleIDs += ubbFlagID;
   obstacleIDs += ubbInflowFlagID;
   obstacleIDs += pressureFlagID;
   obstacleIDs += pressureOutflowFlagID;
   obstacleIDs += freeSlipFlagID;
   obstacleIDs += outletFlagID;
   obstacleIDs += movingObstacleFlagID;

   // initialize outflowIDs
   Set< FlagUID > outflowIDs;
   outflowIDs += pressureOutflowFlagID;
   outflowIDs += outletFlagID;

   // initialize outflowIDs
   Set< FlagUID > inflowIDs;
   inflowIDs += ubbInflowFlagID;

   // initialize freeSlipIDs
   Set< FlagUID > freeSlipIDs;
   freeSlipIDs += freeSlipFlagID;

   // create callable function wrapper with input arguments 1 and 2 unset, whereas arguments 3, 4 and 5 are set to be
   // obstacleIDs, outflowIDs, and inflowIDs, respectively; this is necessary for field::addFlagFieldToStorage()
   auto ffInitFunc = std::bind(internal::flagFieldInitFunction< FlagField_T >, std::placeholders::_1,
                               std::placeholders::_2, obstacleIDs, outflowIDs, inflowIDs, freeSlipIDs);

   // IMPORTANT REMARK: The flag field needs two ghost layers because of function advectMass(). There, the flags of all
   // D3Q* neighbors are determined for each cell, including cells in the first ghost layer. Therefore, all D3Q*
   // neighbors of the first ghost layer must be accessible. This requires a second ghost layer.
   flagFieldID_ = field::addFlagFieldToStorage< FlagField_T >(blockForest, "Flags", uint_c(2), true, ffInitFunc);

   // create FlagInfo
   flagInfo_ = walberla::free_surface::FlagInfo< FlagField_T >(obstacleIDs, outflowIDs, inflowIDs, freeSlipIDs);
   WALBERLA_ASSERT(flagInfo_.isConsistentAcrossBlocksAndProcesses(blockForest, flagFieldID_));

   // add boundary handling to blockForest
   handlingID_ = blockForest_->addBlockData(
      std::make_shared<
         internal::BoundaryBlockDataHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T > >(this),
      "Boundary Handling");

   // create communication object with fill level field, since fill levels determine the flags during the simulation
   comm_.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(fillFieldID_));
}

// define IDs (static const variables)
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::noSlipFlagID =
      field::FlagUID("NoSlip");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::ubbFlagID =
      field::FlagUID("UBB");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::ubbInflowFlagID =
      field::FlagUID("UBB_Inflow");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::pressureFlagID =
      field::FlagUID("Pressure");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::pressureOutflowFlagID =
      field::FlagUID("PressureOutflow");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::outletFlagID =
      field::FlagUID("Outlet");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::freeSlipFlagID =
      field::FlagUID("FreeSlip");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const field::FlagUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::movingObstacleFlagID =
      field::FlagUID("MovingObstacle");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::noSlipBoundaryID =
      BoundaryUID("NoSlip");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::ubbBoundaryID =
      BoundaryUID("UBB");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::ubbInflowBoundaryID =
      BoundaryUID("UBB_Inflow");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::pressureBoundaryID =
      BoundaryUID("Pressure");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                             ParticleAccessor_T >::pressureOutflowBoundaryID =
   BoundaryUID("PressureOutflow");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::outletBoundaryID =
      BoundaryUID("Outlet");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::freeSlipBoundaryID =
      BoundaryUID("FreeSlip");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
const BoundaryUID AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                             ParticleAccessor_T >::movingObstacleBoundaryID =
   BoundaryUID("MovingObstacle");

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
geometry::initializer::InitializationManager
   AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::getInitManager()
{
   using namespace geometry::initializer;

   InitializationManager initManager(blockForest_->getBlockStorage());

   // define initializers
   auto cellIntvInit = std::make_shared< BoundaryFromCellInterval< BoundaryHandling_T > >(*blockForest_, handlingID_);
   auto borderInit   = std::make_shared< BoundaryFromDomainBorder< BoundaryHandling_T > >(*blockForest_, handlingID_);
   auto imgInit =
      std::make_shared< BoundaryFromImage< BoundaryHandling_T, geometry::GrayScaleImage > >(*blockForest_, handlingID_);
   auto bodyInit = std::make_shared< OverlapFieldFromBody >(*blockForest_, fillFieldID_);

   // register initializers
   initManager.registerInitializer("CellInterval", cellIntvInit);
   initManager.registerInitializer("Border", borderInit);
   initManager.registerInitializer("Image", imgInit);
   initManager.registerInitializer("Body", bodyInit);

   return initManager;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::initFromConfig(
   const Config::BlockHandle& configBlock)
{
   // initialize from config file
   getInitManager().init(configBlock);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
template< typename Body_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::addFreeSurfaceObject(
   const Body_T& body, bool addOrSubtract)
{
   geometry::initializer::OverlapFieldFromBody(*blockForest_, fillFieldID_).init(body, addOrSubtract);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setNoSlipAtBorder(
   stencil::Direction d, cell_idx_t wallDistance)
{
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > init(*blockForest_, handlingID_);
   init.init(noSlipFlagID, d, wallDistance);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setNoSlipAtAllBorders(
   cell_idx_t wallDistance)
{
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > init(*blockForest_, handlingID_);
   init.initAllBorders(noSlipFlagID, wallDistance);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setNoSlipInCell(
   const Cell& globalCell)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);

      // transform cell in global coordinates to cell in block local coordinates
      Cell blockLocalCell;
      blockForest_->transformGlobalToBlockLocalCell(blockLocalCell, *blockIt, globalCell);

      handling->forceBoundary(noSlipFlagID, blockLocalCell[0], blockLocalCell[1], blockLocalCell[2]);
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setFreeSlipAtBorder(
   stencil::Direction d, cell_idx_t wallDistance)
{
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > init(*blockForest_, handlingID_);
   init.init(freeSlipFlagID, d, wallDistance);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                ParticleAccessor_T >::setFreeSlipAtAllBorders(cell_idx_t wallDistance)
{
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > init(*blockForest_, handlingID_);
   init.initAllBorders(freeSlipFlagID, wallDistance);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setFreeSlipInCell(
   const Cell& globalCell)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);

      // transform cell in global coordinates to cell in block local coordinates
      Cell blockLocalCell;
      blockForest_->transformGlobalToBlockLocalCell(blockLocalCell, *blockIt, globalCell);

      handling->forceBoundary(freeSlipFlagID, blockLocalCell[0], blockLocalCell[1], blockLocalCell[2]);
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setPressure(
   real_t density)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);
      Pressure_T& pressure =
         handling->template getBoundaryCondition< Pressure_T >(handling->getBoundaryUID(pressureFlagID));
      pressure.setLatticeDensity(density);
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setUBBInCell(
   const Cell& globalCell, const Vector3< real_t >& velocity)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);

      typename UBB_Inflow_T::Velocity ubbVel(velocity);

      // transform cell in global coordinates to cell in block-local coordinates
      Cell blockLocalCell;
      blockForest_->transformGlobalToBlockLocalCell(blockLocalCell, *blockIt, globalCell);

      // get block cell bounding box to check if cell is contained in block
      CellInterval blockCellBB = blockForest_->getBlockCellBB(*blockIt);

      // flag field has two ghost layers so blockCellBB is actually larger than returned; this is relevant for setups
      // where the UBB is set in a ghost layer cell
      blockCellBB.expand(cell_idx_c(2));

      if (blockCellBB.contains(globalCell))
      {
         handling->forceBoundary(ubbFlagID, blockLocalCell[0], blockLocalCell[1], blockLocalCell[2], ubbVel);
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setInflowInCell(
   const Cell& globalCell, const Vector3< real_t >& velocity)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);

      typename UBB_Inflow_T::Velocity ubbVel(velocity);

      // transform cell in global coordinates to cell in block-local coordinates
      Cell blockLocalCell;
      blockForest_->transformGlobalToBlockLocalCell(blockLocalCell, *blockIt, globalCell);

      // get block cell bounding box to check if cell is contained in block
      CellInterval blockCellBB = blockForest_->getBlockCellBB(*blockIt);

      // flag field has two ghost layers so blockCellBB is actually larger than returned; this is relevant for setups
      // where the UBB is set in a ghost layer cell
      blockCellBB.expand(cell_idx_c(2));

      if (blockCellBB.contains(globalCell))
      {
         handling->forceBoundary(ubbInflowFlagID, blockLocalCell[0], blockLocalCell[1], blockLocalCell[2], ubbVel);
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::setPressureOutflow(
   real_t density)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);
      Pressure_T& pressure =
         handling->template getBoundaryCondition< Pressure_T >(handling->getBoundaryUID(pressureOutflowFlagID));
      pressure.setLatticeDensity(density);
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >::enableBubbleOutflow(
   BubbleModelBase* bubbleModel)
{
   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      BoundaryHandling_T* const handling = blockIt->template getData< BoundaryHandling_T >(handlingID_);

      // get pressure from boundary handling
      Pressure_T& pressure =
         handling->template getBoundaryCondition< Pressure_T >(handling->getBoundaryUID(pressureFlagID));
      Pressure_T& pressureOutflow =
         handling->template getBoundaryCondition< Pressure_T >(handling->getBoundaryUID(pressureOutflowFlagID));

      // set pressure in bubble model
      pressure.setBubbleModel(bubbleModel);
      pressureOutflow.setBubbleModel(bubbleModel);
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
Vector3< bool > AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                           ParticleAccessor_T >::isObstacleInGlobalGhostLayer()
{
   Vector3< bool > isObstacleInGlobalGhostLayer(false, false, false);

   for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
   {
      const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID_);

      const CellInterval domainCellBB = blockForest_->getDomainCellBB();

      // disable OpenMP such that loop termination works correctly
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, uint_c(1), omp critical, {
         // get cell in global coordinates
         Cell globalCell = Cell(x, y, z);
         blockForest_->transformBlockLocalToGlobalCell(globalCell, *blockIt);

         // check if the current cell is located in a global ghost layer
         const bool isCellInGlobalGhostLayerX =
            globalCell[0] < domainCellBB.xMin() || globalCell[0] > domainCellBB.xMax();

         const bool isCellInGlobalGhostLayerY =
            globalCell[1] < domainCellBB.yMin() || globalCell[1] > domainCellBB.yMax();

         const bool isCellInGlobalGhostLayerZ =
            globalCell[2] < domainCellBB.zMin() || globalCell[2] > domainCellBB.zMax();

         // skip corners, as they do not influence periodic communication
         if ((isCellInGlobalGhostLayerX && (isCellInGlobalGhostLayerY || isCellInGlobalGhostLayerZ)) ||
             (isCellInGlobalGhostLayerY && isCellInGlobalGhostLayerZ))
         {
            continue;
         }

         if (!isObstacleInGlobalGhostLayer[0] && isCellInGlobalGhostLayerX &&
             isPartOfMaskSet(flagField->get(x, y, z), flagField->getMask(flagInfo_.getObstacleIDSet())))
         {
            isObstacleInGlobalGhostLayer[0] = true;
         }

         if (!isObstacleInGlobalGhostLayer[1] && isCellInGlobalGhostLayerY &&
             isPartOfMaskSet(flagField->get(x, y, z), flagField->getMask(flagInfo_.getObstacleIDSet())))
         {
            isObstacleInGlobalGhostLayer[1] = true;
         }

         if (!isObstacleInGlobalGhostLayer[2] && isCellInGlobalGhostLayerZ &&
             isPartOfMaskSet(flagField->get(x, y, z), flagField->getMask(flagInfo_.getObstacleIDSet())))
         {
            isObstacleInGlobalGhostLayer[2] = true;
         }

         if (isObstacleInGlobalGhostLayer[0] && isObstacleInGlobalGhostLayer[1] && isObstacleInGlobalGhostLayer[2])
         {
            break; // there is no need to check other cells on this block
         }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP
   }

   mpi::allReduceInplace(isObstacleInGlobalGhostLayer, mpi::LOGICAL_OR);

   return isObstacleInGlobalGhostLayer;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename ParticleAccessor_T >
void AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T,
                                ParticleAccessor_T >::initFlagsFromFillLevel()
{
   const Vector3< bool > isObstacleInGlobalGhostLayer = this->isObstacleInGlobalGhostLayer();

   WALBERLA_ROOT_SECTION()
   {
      if ((blockForest_->isXPeriodic() && isObstacleInGlobalGhostLayer[0]) ||
          (blockForest_->isYPeriodic() && isObstacleInGlobalGhostLayer[1]) ||
          (blockForest_->isZPeriodic() && isObstacleInGlobalGhostLayer[2]))
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
            "WARNING: An obstacle cell is located in a global outermost ghost layer in a periodic "
            "direction. Be aware that due to periodicity, this obstacle cell will be "
            "overwritten during communication.");
      }
   }

   // communicate fill level (neighborhood is used in initialization)
   comm_();

   // initialize fill level in boundaries (with value 1), i.e., obstacles such that the bubble model does not detect
   // obstacles as gas cells
   walberla::free_surface::initFillLevelsInBoundaries< BoundaryHandling_T, typename LatticeModel_T::Stencil,
                                                       ScalarField_T >(blockForest_, handlingID_, fillFieldID_);

   // clear and initialize flags in every cell according to the fill level
   walberla::free_surface::initFlagsFromFillLevels< BoundaryHandling_T, typename LatticeModel_T::Stencil, FlagField_T,
                                                    const ScalarField_T >(blockForest_, flagInfo_, handlingID_,
                                                                          fillFieldID_);
}

} // namespace free_surface
} // namespace antidunes
} // namespace walberla
