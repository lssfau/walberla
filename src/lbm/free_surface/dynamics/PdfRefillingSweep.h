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
//! \file PdfRefillingSweep.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \author Michael Zikeli
//! \brief Sweeps for refilling cells (i.e. reinitializing PDFs) after the cell was converted from gas to interface.
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/StringUtility.h"
#include "core/cell/Cell.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/IBlock.h"

#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/dynamics/PdfRefillingModel.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"

#include "stencil/D3Q6.h"

#include <functional>
#include <vector>

#include "PdfRefillingModel.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Base class for all sweeps to reinitialize (refill) all PDFs in cells that were converted from gas to interface.
 * This is required since gas cells do not have PDFs. The sweep expects that a previous sweep has set the
 * "convertedFromGasToInterface" flag and reinitializes the PDFs in all cells with this flag according to the specified
 * model.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T >
class RefillingSweepBase
{
 public:
   using PdfField_T = lbm::PdfField< LatticeModel_T >;
   using flag_t     = typename FlagField_T::flag_t;
   using Stencil_T  = typename LatticeModel_T::Stencil;

   RefillingSweepBase(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                      const FlagInfo< FlagField_T >& flagInfo, bool useDataFromGhostLayers)
      : pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID), flagInfo_(flagInfo),
        useDataFromGhostLayers_(useDataFromGhostLayers)
   {}

   virtual void operator()(IBlock* const block) = 0;

   virtual ~RefillingSweepBase() = default;

   real_t getAverageDensityAndVelocity(const Cell& cell, const PdfField_T& pdfField, const FlagField_T& flagField,
                                       const FlagInfo< FlagField_T >& flagInfo, Vector3< real_t >& avgVelocity)
   {
      std::vector< bool > validStencilIndices(Stencil_T::Size, false);
      return getAverageDensityAndVelocity(cell, pdfField, flagField, flagInfo, avgVelocity, validStencilIndices);
   }

   // also stores stencil indices of valid neighboring cells (liquid/ non-newly converted interface) in a vector
   real_t getAverageDensityAndVelocity(const Cell& cell, const PdfField_T& pdfField, const FlagField_T& flagField,
                                       const FlagInfo< FlagField_T >& flagInfo, Vector3< real_t >& avgVelocity,
                                       std::vector< bool >& validStencilIndices);

   // returns the averaged PDFs of valid neighboring cells (liquid/ non-newly converted interface)
   std::vector< real_t > getAveragePdfs(const Cell& cell, const PdfField_T& pdfField, const FlagField_T& flagField,
                                        const FlagInfo< FlagField_T >& flagInfo);

 protected:
   BlockDataID pdfFieldID_;
   ConstBlockDataID flagFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
   bool useDataFromGhostLayers_;
}; // class RefillingSweepBase

/***********************************************************************************************************************
 * Base class for refilling models that need to obtain information by extrapolation from neighboring cells.
 *
 * The parameter "useDataFromGhostLayers" is useful, when reducedCommunication is used, i.e., when not all PDFs are
 * communicated in the PDF field but only those that are actually required in a certain direction. This is currently not
 * used in the free surface part of waLBerla: In SurfaceDynamicsHandler, SimpleCommunication uses the default PackInfo
 * of the GhostLayerField for PDF communication. The optimized version would be using the PdfFieldPackInfo (in
 * src/lbm/communication).
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExtrapolationRefillingSweepBase : public RefillingSweepBase< LatticeModel_T, FlagField_T >
{
 public:
   using RefillingSweepBase_T = RefillingSweepBase< LatticeModel_T, FlagField_T >;
   using PdfField_T           = typename RefillingSweepBase_T::PdfField_T;
   using flag_t               = typename RefillingSweepBase_T::flag_t;
   using Stencil_T            = typename RefillingSweepBase_T::Stencil_T;

   ExtrapolationRefillingSweepBase(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                                   const ConstBlockDataID& fillFieldID, const FlagInfo< FlagField_T >& flagInfo,
                                   uint_t extrapolationOrder, bool useDataFromGhostLayers)
      : RefillingSweepBase_T(pdfFieldID, flagFieldID, flagInfo, useDataFromGhostLayers), fillFieldID_(fillFieldID),
        extrapolationOrder_(extrapolationOrder)
   {}

   virtual ~ExtrapolationRefillingSweepBase() = default;

   virtual void operator()(IBlock* const block) = 0;

   /********************************************************************************************************************
    * Find the lattice direction in the given stencil that corresponds best to the provided direction.
    *
    * Mostly copied from src/lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h.
    *******************************************************************************************************************/
   Vector3< cell_idx_t > findCorrespondingLatticeDirection(const Vector3< real_t >& direction);

   /********************************************************************************************************************
    * The extrapolation direction is chosen such that it most closely resembles the surface normal in the cell.
    *
    * The normal has to be recomputed here and MUST NOT be taken directly from the normal field because the fill levels
    * will have changed since the last computation of the normal. As the normal is computed from the fill levels, it
    * must be recomputed to have an up-to-date normal.
    *******************************************************************************************************************/
   Vector3< cell_idx_t > findExtrapolationDirection(const Cell& cell, const FlagField_T& flagField,
                                                    const ScalarField_T& fillField);

   /********************************************************************************************************************
    * Determine the number of applicable (liquid or interface) cells for extrapolation in the given extrapolation
    * direction.
    *******************************************************************************************************************/
   uint_t getNumberOfExtrapolationCells(const Cell& cell, const FlagField_T& flagField, const PdfField_T& pdfField,
                                        const Vector3< cell_idx_t >& extrapolationDirection);

   /********************************************************************************************************************
    * Get the non-equilibrium part of all PDFs in "cell" and store them in a std::vector<real_t>.
    *******************************************************************************************************************/
   std::vector< real_t > getNonEquilibriumPdfsInCell(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField);

   /********************************************************************************************************************
    * Get all PDFs in "cell" and store them in a std::vector<real_t>.
    *******************************************************************************************************************/
   std::vector< real_t > getPdfsInCell(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField);

   /********************************************************************************************************************
    * Set the PDFs in cell "x" according to the following linear combination:
    *   f(x,q) = f(x,q) + 3 * f^{get}(x+e,q) - 3 * f^{get}(x+2e,q) + 1 * f^{get}(x+3e,q)
    *       with x: cell position
    *            q: index of the respective PDF
    *            e: extrapolation direction
    *            f^{get}: PDF specified by getPdfFunc
    *       "includeThisCell" defines whether f(x,q) is included
    *
    * Note: this function does NOT assert out of bounds access, i.e., it may only be called for cells that have at least
    *       three neighboring cells in extrapolationDirection.
    *******************************************************************************************************************/
   void applyQuadraticExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc);

   /********************************************************************************************************************
    * Set the PDFs in cell (x) according to the following linear combination:
    *   f(x,q) = f(x,q) + 2 * f^{get}(x+1e,q) - 1 * f^{get}(x+2e,q)
    *       with x: cell position
    *            q: index of the respective PDF
    *            e: extrapolation direction
    *            f^{get}: PDF specified by getPdfFunc
    *       "includeThisCell" defines whether f(x,q) is included
    *
    * Note: this function does NOT assert out of bounds access, i.e., it may only be called for cells that have at least
    *       two neighboring cells in extrapolationDirection.
    *******************************************************************************************************************/
   void applyLinearExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc);

   /********************************************************************************************************************
    * Set the PDFs in cell (x) according to the following linear combination:
    *   f(x,q) = f(x,q) + 1 * f^{get}(x+1e,q)
    *       with x: cell position
    *            q: index of the respective PDF
    *            e: extrapolation direction
    *            f^{get}: PDF specified by getPdfFunc
    *       "includeThisCell" defines whether f(x,q) is included
    *
    * Note: this function does NOT assert out of bounds access, i.e., it may only be called for cells that have at least
    *       a neighboring cell in extrapolationDirection.
    *******************************************************************************************************************/
   void applyConstantExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc);

 protected:
   ConstBlockDataID fillFieldID_;
   uint_t extrapolationOrder_;
}; // class ExtrapolationRefillingSweepBase

/***********************************************************************************************************************
 * EquilibriumRefillingSweep:
 * PDFs are initialized with the equilibrium based on the average density and velocity from neighboring liquid and
 * (non-newly converted) interface cells.
 *
 * Reference: dissertation of N. Thuerey, 2007, section 4.3
 *
 *  f(x,q) = f^{eq}(x+e,q)
 *      with x: cell position
 *           q: index of the respective PDF
 *           e: direction of a valid neighbor
 **********************************************************************************************************************/

template< typename LatticeModel_T, typename FlagField_T >
class EquilibriumRefillingSweep : public RefillingSweepBase< LatticeModel_T, FlagField_T >
{
 public:
   using RefillingSweepBase_T = RefillingSweepBase< LatticeModel_T, FlagField_T >;
   using PdfField_T           = typename RefillingSweepBase_T::PdfField_T;
   using flag_t               = typename RefillingSweepBase_T::flag_t;
   using Stencil_T            = typename RefillingSweepBase_T::Stencil_T;

   EquilibriumRefillingSweep(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                             const FlagInfo< FlagField_T >& flagInfo, bool useDataFromGhostLayers)
      : RefillingSweepBase_T(pdfFieldID, flagFieldID, flagInfo, useDataFromGhostLayers)
   {}

   ~EquilibriumRefillingSweep() override = default;

   void operator()(IBlock* const block) override;
}; // class EquilibriumRefillingSweep

/***********************************************************************************************************************
 * AverageRefillingSweep:
 * PDFs are initialized with the average of the PDFs in the same direction from applicable neighboring cells.
 *
 *  f_i(x,q) = \sum_N( f_i(x+e,q) ) / N
 *      with x: cell position
 *           i: PDF index, i.e., direction
 *           q: index of the respective PDF
 *           e: direction of a valid neighbor
 *           N: number of applicable neighbors
 *           sum_N: sum over all applicable neighbors
 *
 * Reference: not available in literature (as of 06/2022).
 **********************************************************************************************************************/

template< typename LatticeModel_T, typename FlagField_T >
class AverageRefillingSweep : public RefillingSweepBase< LatticeModel_T, FlagField_T >
{
 public:
   using RefillingSweepBase_T = RefillingSweepBase< LatticeModel_T, FlagField_T >;
   using PdfField_T           = typename RefillingSweepBase_T::PdfField_T;
   using flag_t               = typename RefillingSweepBase_T::flag_t;
   using Stencil_T            = typename RefillingSweepBase_T::Stencil_T;

   AverageRefillingSweep(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                         const FlagInfo< FlagField_T >& flagInfo, bool useDataFromGhostLayers)
      : RefillingSweepBase_T(pdfFieldID, flagFieldID, flagInfo, useDataFromGhostLayers)
   {}

   ~AverageRefillingSweep() override = default;

   void operator()(IBlock* const block) override;
}; // class AverageRefillingSweep

/***********************************************************************************************************************
 * EquilibriumAndNonEquilibriumRefilling:
 * First reconstruct the equilibrium values according to the "EquilibriumRefilling". Then extrapolate the
 * non-equilibrium part of the PDFs in the direction of the surface normal and add it to the (equilibrium-)
 * reinitialized PDFs.
 *
 * Reference: equation (51) in  Peng et al., "Implementation issues and benchmarking of lattice Boltzmann method for
 *            moving rigid particle simulations in a viscous flow", 2015, doi: 10.1016/j.camwa.2015.08.027)
 *
 *  f_q(x,t+dt) = f_q^{eq}(x,t) + f_q^{neq}(x+e*dt,t+dt)
 *      with x: cell position
 *           q: index of the respective PDF
 *           e: direction of a valid neighbor/ extrapolation direction
 *           t: current time step
 *           dt: time step width
 *           f^{eq}: equilibrium PDF with velocity and density averaged from neighboring cells
 *           f^{neq}: non-equilibrium PDF (f - f^{eq})
 *
 * Note: Analogously as in the "ExtrapolationRefilling", the expression "f_q^{neq}(x+e*dt,t+dt)" can also be obtained by
 *       extrapolation (in literature, only zeroth order extrapolation is used/documented):
 *          - zeroth order: f_q^{neq}(x+e*dt,t+dt)
 *          - first order: 2 * f_q^{neq}(x+e*dt,t+dt) - 1 * f_q^{neq}(x+2e*dt,t+dt)
 *          - second order: 3 * f_q^{neq}(x+e*dt,t+dt) - 3 * f_q^{neq}(x+2e*dt,t+dt) + f_q^{neq}(x+3e*dt,t+dt)
 * If not enough cells are available for the chosen extrapolation order, the algorithm falls back to the corresponding
 * lower order. If even zeroth order can not be applied, only f_q^{eq}(x,t) is considered which corresponds to
 * "EquilibriumRefilling".
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class EquilibriumAndNonEquilibriumRefillingSweep
   : public ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
{
 public:
   using ExtrapolationRefillingSweepBase_T =
      ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;
   using RefillingSweepBase_T = typename ExtrapolationRefillingSweepBase_T::RefillingSweepBase_T;
   using PdfField_T           = typename ExtrapolationRefillingSweepBase_T::PdfField_T;
   using flag_t               = typename ExtrapolationRefillingSweepBase_T::flag_t;
   using Stencil_T            = typename ExtrapolationRefillingSweepBase_T::Stencil_T;

   EquilibriumAndNonEquilibriumRefillingSweep(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                                              const ConstBlockDataID& fillFieldID,
                                              const FlagInfo< FlagField_T >& flagInfo, uint_t extrapolationOrder,
                                              bool useDataFromGhostLayers)
      : ExtrapolationRefillingSweepBase_T(pdfFieldID, flagFieldID, fillFieldID, flagInfo, extrapolationOrder,
                                          useDataFromGhostLayers)
   {}

   ~EquilibriumAndNonEquilibriumRefillingSweep() override = default;

   void operator()(IBlock* const block) override;
}; // class EquilibriumAndNonEquilibriumRefillingSweep

/***********************************************************************************************************************
 * ExtrapolationRefilling:
 * Extrapolate the PDFs of one or more cells in the direction of the surface normal.
 *
 * Reference: equation (50) in  Peng et al., "Implementation issues and benchmarking of lattice Boltzmann method for
 *            moving rigid particle simulations in a viscous flow", 2015, doi: 10.1016/j.camwa.2015.08.027)
 *
 *  f_q(x,t+dt) = 3 * f_q(x+e*dt,t+dt) - 3 * f_q(x+2e*dt,t+dt) + f_q(x+3e*dt,t+dt)
 *        with x: cell position
 *             q: index of the respective PDF
 *             e: direction of a valid neighbor/ extrapolation direction
 *             t: current time step
 *             dt: time step width
 * Note: The equation contains a second order extrapolation, however other options are also available. If not enough
 *       cells are available for second order extrapolation, the algorithm falls back to the next applicable
 *       lower order:
 *           - second order: 3 * f_q(x+e*dt,t+dt) - 3 * f_q(x+2e*dt,t+dt) + f_q(x+3e*dt,t+dt)
 *           - first order: 2 * f_q(x+e*dt,t+dt) - 1 * f_q(x+2e*dt,t+dt)
 *           - zeroth order: f_q(x+e*dt,t+dt)
 * If even zeroth order can not be applied, only f_q^{eq}(x,t) is considered which corresponds to
 * "EquilibriumRefilling".
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class ExtrapolationRefillingSweep
   : public ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
{
 public:
   using ExtrapolationRefillingSweepBase_T =
      ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;
   using RefillingSweepBase_T = typename ExtrapolationRefillingSweepBase_T::RefillingSweepBase_T;
   using PdfField_T           = typename ExtrapolationRefillingSweepBase_T::PdfField_T;
   using flag_t               = typename ExtrapolationRefillingSweepBase_T::flag_t;
   using Stencil_T            = typename ExtrapolationRefillingSweepBase_T::Stencil_T;

   ExtrapolationRefillingSweep(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                               const ConstBlockDataID& fillFieldID, const FlagInfo< FlagField_T >& flagInfo,
                               uint_t extrapolationOrder, bool useDataFromGhostLayers)
      : ExtrapolationRefillingSweepBase_T(pdfFieldID, flagFieldID, fillFieldID, flagInfo, extrapolationOrder,
                                          useDataFromGhostLayers)
   {}

   ~ExtrapolationRefillingSweep() override = default;

   void operator()(IBlock* const block) override;
}; // class ExtrapolationRefillingSweep

/***********************************************************************************************************************
 * GradsMomentsRefilling:
 * Reconstruct missing PDFs based on Grad's moment closure.
 *
 * References: - equation (11) in Chikatamarla et al., "Grad’s approximation for missing data in lattice Boltzmann
 *               simulations", 2006, doi: 10.1209/epl/i2005-10535-x
 *             - equation (10) in Dorscher et al., "Grad’s approximation for moving and stationary walls in entropic
 *               lattice Boltzmann simulations", 2015, doi: 10.1016/j.jcp.2015.04.017
 *
 * The following equation is a rewritten and easier version of the equation in the above references:
 *  f_q(x,t+dt) = f_q^{eq}(x,t) +
 *                w_q * rho / 2 / cs^2 / omega * (du_a / dx_b + du_b / dx_a)(cs^2 * delta_{ab} - c_{q,a}c_{q,b} )
 *      with x: cell position
 *           q: index of the respective PDF
 *           t: current time step
 *           f^{eq}: equilibrium PDF with velocity and density averaged from neighboring cells
 *           w_q: lattice weight
 *           rho: density averaged from neighboring cells
 *           cs: lattice speed of sound
 *           omega: relaxation rate
 *           du_a / dx_b: gradient of the velocity (in index notation)
 *           delta_{ab}: Kronecker delta (in index notation)
 *           c_q{q,a}: lattice velocity (in index notation)
 *
 * The velocity gradient is computed using a first order central finite difference scheme if two neighboring cells
 * are available. Otherwise, a first order upwind scheme is applied.
 * IMPORTANT REMARK: The current implementation only works for dx=1 (this is assumed in the gradient computation).
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T >
class GradsMomentsRefillingSweep : public RefillingSweepBase< LatticeModel_T, FlagField_T >
{
 public:
   using RefillingSweepBase_T = RefillingSweepBase< LatticeModel_T, FlagField_T >;
   using PdfField_T           = typename RefillingSweepBase_T::PdfField_T;
   using flag_t               = typename RefillingSweepBase_T::flag_t;
   using Stencil_T            = typename RefillingSweepBase_T::Stencil_T;

   GradsMomentsRefillingSweep(const BlockDataID& pdfFieldID, const ConstBlockDataID& flagFieldID,
                              const FlagInfo< FlagField_T >& flagInfo, real_t relaxRate, bool useDataFromGhostLayers)

      : RefillingSweepBase_T(pdfFieldID, flagFieldID, flagInfo, useDataFromGhostLayers), relaxRate_(relaxRate)
   {}

   ~GradsMomentsRefillingSweep() override = default;

   void operator()(IBlock* const block) override;

   // compute the gradient of the velocity in the specified direction
   // - using a first order central finite difference scheme if two valid neighboring cells are available
   // - using a first order upwind scheme if only one valid neighboring cell is available
   // - assuming a gradient of zero if no valid neighboring cell is available
   Vector3< real_t > getVelocityGradient(stencil::Direction direction, const Cell& cell, const PdfField_T* pdfField,
                                         const Vector3< real_t >& avgVelocity,
                                         const std::vector< bool >& validStencilIndices);

 private:
   real_t relaxRate_;
}; // class GradsMomentApproximationRefilling

} // namespace free_surface
} // namespace walberla

#include "PdfRefillingSweep.impl.h"