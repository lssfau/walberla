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
//! \file Boundary.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Boundary handling for the free surface LBM module.
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "field/GhostLayerField.h"

#include "geometry/initializer/InitializationManager.h"
#include "geometry/initializer/OverlapFieldFromBody.h"

#include "lbm/boundary/all.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/InitFunctions.h"
#include "lbm/free_surface/InterfaceFromFillLevel.h"
#include "lbm/free_surface/boundary/SimplePressureWithFreeSurface.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Boundary handling for the free surface LBM extension.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T >
class FreeSurfaceBoundaryHandling
{
 public:
   using flag_t    = typename FlagField_T::value_type;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using CommunicationStencil_T =
      typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   // boundary
   using NoSlip_T   = lbm::NoSlip< LatticeModel_T, flag_t >;
   using FreeSlip_T = lbm::FreeSlip< LatticeModel_T, FlagField_T >;
   using UBB_T      = lbm::UBB< LatticeModel_T, flag_t >;
   using Pressure_T = SimplePressureWithFreeSurface< LatticeModel_T, FlagField_T >;
   using Outlet_T   = lbm::Outlet< LatticeModel_T, FlagField_T, 4, 3 >;
   using UBB_Inflow_T =
      lbm::UBB< LatticeModel_T, flag_t >; // creates interface cells in the direction of the prescribed velocity, i.e.,
                                          // represents an inflow boundary condition

   // handling type
   using BoundaryHandling_T =
      BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, UBB_T, UBB_Inflow_T, Pressure_T, Pressure_T, Outlet_T,
                        FreeSlip_T >; // 2 pressure boundaries with different densities, e.g., inflow-outflow

   using FlagInfo_T = FlagInfo< FlagField_T >;

   // constructor
   FreeSurfaceBoundaryHandling(const std::shared_ptr< StructuredBlockForest >& blockForest, BlockDataID pdfFieldID,
                               BlockDataID fillLevelID);

   // initialize fluid field from config file using (quotes indicate the string to be used in the file):
   // - "CellInterval" (see src/geometry/initializer/BoundaryFromCellInterval.h)
   // - "Border" (see src/geometry/initializer/BoundaryFromDomainBorder.h)
   // - "Image" (see src/geometry/initializer/BoundaryFromImage.h)
   // - "Body" (see src/geometry/initializer/OverlapFieldFromBody.h)
   inline void initFromConfig(const Config::BlockHandle& block);

   // initialize free surface object from geometric body (see src/geometry/initializer/OverlapFieldFromBody.h)
   template< typename Body_T >
   inline void addFreeSurfaceObject(const Body_T& body, bool addOrSubtract = false);

   // clear and initialize flags in every cell according to the fill level and initialize fill level in boundaries (with
   // value 1) and obstacles such that the bubble model does not detect obstacles as gas cells
   void initFlagsFromFillLevel();

   inline void setNoSlipAtBorder(stencil::Direction d, cell_idx_t wallDistance = cell_idx_c(0));
   inline void setNoSlipAtAllBorders(cell_idx_t wallDistance = cell_idx_c(0));
   void setNoSlipInCell(const Cell& globalCell);

   inline void setFreeSlipAtBorder(stencil::Direction d, cell_idx_t wallDistance = cell_idx_c(0));
   inline void setFreeSlipAtAllBorders(cell_idx_t wallDistance = cell_idx_c(0));
   void setFreeSlipInCell(const Cell& globalCell);

   void setUBBInCell(const Cell& globalCell, const Vector3< real_t >& velocity);

   // UBB that generates interface cells to resemble an inflow boundary
   void setInflowInCell(const Cell& globalCell, const Vector3< real_t >& velocity);

   inline void setPressure(real_t density);
   void setPressureOutflow(real_t density);
   void setBodyForce(const Vector3< real_t >& bodyForce);

   void enableBubbleOutflow(BubbleModelBase* bubbleModel);

   // checks if an obstacle cell is located in an outermost ghost layer (corners are explicitly ignored, as they do not
   // influence periodic communication)
   Vector3< bool > isObstacleInGlobalGhostLayer();

   // flag management
   const FlagInfo< FlagField_T >& getFlagInfo() const { return flagInfo_; }

   // flag IDs
   static const field::FlagUID noSlipFlagID;
   static const field::FlagUID ubbFlagID;
   static const field::FlagUID ubbInflowFlagID;
   static const field::FlagUID pressureFlagID;
   static const field::FlagUID pressureOutflowFlagID;
   static const field::FlagUID outletFlagID;
   static const field::FlagUID freeSlipFlagID;

   // boundary IDs
   static const BoundaryUID noSlipBoundaryID;
   static const BoundaryUID ubbBoundaryID;
   static const BoundaryUID ubbInflowBoundaryID;
   static const BoundaryUID pressureBoundaryID;
   static const BoundaryUID pressureOutflowBoundaryID;
   static const BoundaryUID outletBoundaryID;
   static const BoundaryUID freeSlipBoundaryID;

   inline BlockDataID getHandlingID() const { return handlingID_; }
   inline BlockDataID getPdfFieldID() const { return pdfFieldID_; }
   inline BlockDataID getFillFieldID() const { return fillFieldID_; }
   inline BlockDataID getFlagFieldID() const { return flagFieldID_; }

   // executes standard waLBerla boundary handling
   class ExecuteBoundaryHandling
   {
    public:
      ExecuteBoundaryHandling(const BlockDataID& collection) : handlingID_(collection) {}
      void operator()(IBlock* const block) const
      {
         BoundaryHandling_T* const handling = block->getData< BoundaryHandling_T >(handlingID_);
         // reset "near boundary" flags
         handling->refresh();
         (*handling)();
      }

    protected:
      BlockDataID handlingID_;
   }; // class ExecuteBoundaryHandling

   ExecuteBoundaryHandling getBoundarySweep() const { return ExecuteBoundaryHandling(getHandlingID()); }

 private:
   FlagInfo< FlagField_T > flagInfo_;

   // register standard waLBerla initializers
   geometry::initializer::InitializationManager getInitManager();

   std::shared_ptr< StructuredBlockForest > blockForest_;

   BlockDataID flagFieldID_;
   BlockDataID pdfFieldID_;
   BlockDataID fillFieldID_;

   BlockDataID handlingID_;

   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > comm_;
}; // class FreeSurfaceBoundaryHandling

} // namespace free_surface
} // namespace walberla

#include "FreeSurfaceBoundaryHandling.impl.h"
