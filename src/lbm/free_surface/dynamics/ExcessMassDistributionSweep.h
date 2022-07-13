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
//! \file ExcessMassDistributionSweep.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Distribute excess mass, i.e., mass that is undistributed after conversions from interface to liquid or gas.
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FieldClone.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"

#include "ExcessMassDistributionModel.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Distribute excess mass, i.e., mass that is undistributed after cells have been converted from interface to
 * gas/liquid. For example, when converting an interface cell with fill level 1.1 to liquid with fill level 1.0, an
 * excessive mass corresponding to the fill level 0.1 must be distributed to conserve mass.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExcessMassDistributionSweepBase
{
 public:
   ExcessMassDistributionSweepBase(const ExcessMassDistributionModel& excessMassDistributionModel,
                                   BlockDataID fillFieldID, ConstBlockDataID flagFieldID, ConstBlockDataID pdfFieldID,
                                   const FlagInfo< FlagField_T >& flagInfo)
      : excessMassDistributionModel_(excessMassDistributionModel), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        pdfFieldID_(pdfFieldID), flagInfo_(flagInfo)
   {}

   virtual void operator()(IBlock* const block) = 0;

   virtual ~ExcessMassDistributionSweepBase() = default;

 protected:
   /********************************************************************************************************************
    * Determines the number of a cell's
    *  - neighboring newly-converted interface cells
    *  - neighboring interface cells (regardless if newly converted or not)
    *******************************************************************************************************************/
   void getNumberOfInterfaceNeighbors(const FlagField_T* flagField, const Cell& cell, uint_t& newInterfaceNeighbors,
                                      uint_t& interfaceNeighbors);

   /********************************************************************************************************************
    * Determines the number of a cell's neighboring liquid and interface cells.
    *******************************************************************************************************************/
   void getNumberOfEvenlyLiquidAndAllInterfacePreferInterfaceNeighbors(const FlagField_T* flagField, const Cell& cell,
                                                                       uint_t& liquidNeighbors,
                                                                       uint_t& interfaceNeighbors);

   ExcessMassDistributionModel excessMassDistributionModel_;
   BlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID pdfFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
}; // class ExcessMassDistributionSweep

/***********************************************************************************************************************
 * Distribute the excess mass evenly among either
 *  - all neighboring interface cells (see dissertations of T. Pohl, S. Donath, S. Bogner).
 *  - newly converted interface cells (see Koerner et al., 2005)
 *  - old, i.e., non-newly converted interface cells
 *
 * If either no newly converted interface cell or old interface cell is available in the neighborhood, the
 * respective other approach is used as fallback.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExcessMassDistributionSweepInterfaceEvenly
   : public ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
{
 public:
   using ExcessMassDistributionSweepBase_T =
      ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

   ExcessMassDistributionSweepInterfaceEvenly(const ExcessMassDistributionModel& excessMassDistributionModel,
                                              BlockDataID fillFieldID, ConstBlockDataID flagFieldID,
                                              ConstBlockDataID pdfFieldID, const FlagInfo< FlagField_T >& flagInfo)
      : ExcessMassDistributionSweepBase_T(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo)
   {}

   ~ExcessMassDistributionSweepInterfaceEvenly() override = default;

   void operator()(IBlock* const block) override;

 private:
   template< typename PdfField_T >
   void distributeMassEvenly(ScalarField_T* fillField, const FlagField_T* flagField, const PdfField_T* pdfField,
                             const Cell& cell, real_t excessFill);
}; // class ExcessMassDistributionSweepInterfaceEvenly

/***********************************************************************************************************************
 * Distribute the excess mass weighted with the direction of the interface normal among either
 *  - all neighboring interface cells (see section 4.3 in dissertation of N. Thuerey, 2007)
 *  - newly converted interface cells
 *  - old, i.e., non-newly converted interface cells
 *
 * If either no newly converted interface cell or old interface cell is available in the neighborhood, the
 * respective other approach is used as fallback.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExcessMassDistributionSweepInterfaceWeighted
   : public ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
{
 public:
   using ExcessMassDistributionSweepBase_T =
      ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

   ExcessMassDistributionSweepInterfaceWeighted(const ExcessMassDistributionModel& excessMassDistributionModel,
                                                BlockDataID fillFieldID, ConstBlockDataID flagFieldID,
                                                ConstBlockDataID pdfFieldID, const FlagInfo< FlagField_T >& flagInfo,
                                                ConstBlockDataID normalFieldID)
      : ExcessMassDistributionSweepBase_T(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo),
        normalFieldID_(normalFieldID)
   {}

   ~ExcessMassDistributionSweepInterfaceWeighted() override = default;

   void operator()(IBlock* const block) override;

 private:
   template< typename PdfField_T >
   void distributeMassWeighted(ScalarField_T* fillField, const FlagField_T* flagField, const PdfField_T* pdfField,
                               const VectorField_T* normalField, const Cell& cell, bool isNewLiquid, real_t excessFill);

   /********************************************************************************************************************
    * Returns vector with weights for excess mass distribution among neighboring cells.
    *******************************************************************************************************************/
   void getExcessMassWeights(const FlagField_T* flagField, const VectorField_T* normalField, const Cell& cell,
                             bool isNewLiquid, bool useWeightedOld, bool useWeightedAll, bool useWeightedNew,
                             std::vector< real_t >& weights);

   /********************************************************************************************************************
    * Computes the weights for distributing the excess mass based on the direction of the interface normal (see equation
    * (4.9) in dissertation of N. Thuerey, 2007)
    *******************************************************************************************************************/
   void computeWeightWithNormal(real_t n_dot_ci, bool isNewLiquid, typename LatticeModel_T::Stencil::iterator dir,
                                std::vector< real_t >& weights);

   ConstBlockDataID normalFieldID_;

}; // class ExcessMassDistributionSweepInterfaceWeighted

/***********************************************************************************************************************
 * Distribute the excess mass evenly among
 *  - all neighboring liquid and interface cells (see p. 47 in master thesis of M. Lehmann, 2019)
 *  - all neighboring interface cells and only to liquid cells if there exists no neighboring interface cell
 *
 * Neither the fill level, nor the density of liquid cells is modified. Instead, the excess mass is stored in an
 * additional field.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExcessMassDistributionSweepInterfaceAndLiquid
   : public ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
{
 public:
   using ExcessMassDistributionSweepBase_T =
      ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

   ExcessMassDistributionSweepInterfaceAndLiquid(const ExcessMassDistributionModel& excessMassDistributionModel,
                                                 BlockDataID fillFieldID, ConstBlockDataID flagFieldID,
                                                 ConstBlockDataID pdfFieldID, const FlagInfo< FlagField_T >& flagInfo,
                                                 BlockDataID excessMassFieldID)
      : ExcessMassDistributionSweepBase_T(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo),
        excessMassFieldID_(excessMassFieldID), excessMassFieldClone_(excessMassFieldID)
   {}

   ~ExcessMassDistributionSweepInterfaceAndLiquid() override = default;

   void operator()(IBlock* const block) override;

 private:
   template< typename PdfField_T >
   void distributeMassInterfaceAndLiquid(ScalarField_T* fillField, ScalarField_T* dstExcessMassField,
                                         const FlagField_T* flagField, const PdfField_T* pdfField, const Cell& cell,
                                         real_t excessMass);

   BlockDataID excessMassFieldID_;
   field::FieldClone< ScalarField_T, true > excessMassFieldClone_;

}; // class ExcessMassDistributionSweepInterfaceAndLiquid

} // namespace free_surface
} // namespace walberla

#include "ExcessMassDistributionSweep.impl.h"