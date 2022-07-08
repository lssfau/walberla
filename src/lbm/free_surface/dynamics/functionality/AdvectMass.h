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
//! \file AdvectMass.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Calculate the mass advection for a single interface cell.
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingPool.h"

#include "domain_decomposition/IBlock.h"

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"

#include <type_traits>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Calculate the mass advection for a single interface cell and returns the mass delta for this cell.
 *
 * With "useSimpleMassExchange==true", mass exchange is performed according to the paper from Koerner et al., 2005,
 * equation (9). That is, mass is exchanged regardless of the cells' neighborhood.
 *  Mass exchange of interface cell with
 * - obstacle cell: no mass is exchanged
 * - gas cell: no mass is exchanged
 * - liquid cell: mass exchange is determined by the difference between incoming and outgoing PDFs
 * - interface cell: The mass exchange is determined by the difference between incoming and outgoing PDFs weighted with
 *                   the average fill level of both interface cells. This is basically a linear estimate of the wetted
 *                   area between the two cells.
 * - free slip cell: mass is exchanged with the cell where the mirrored PDFs are coming from
 * - inflow cell: mass exchange as with liquid cell
 * - outflow cell: mass exchange as with interface cell (outflow cell is assumed to have the same fill level)
 *
 * With  "useSimpleMassExchange==false", mass is exchanged according to the dissertation of N. Thuerey,
 * 2007, sections 4.1 and 4.4, and table 4.1. However, here, the fill level field and density are used for the
 * computations instead of storing an additional mass field.
 * To ensure a single interface layer, an interface cell is
 * - forced to empty if it has no fluid neighbors.
 * - forced to fill if it has no gas neighbors.
 * This is done by modifying the exchanged PDFs accordingly. If an interface cell is
 * - forced to empty, only outgoing PDFs but no incoming PDFs are allowed.
 * - forced to fill, only incoming PDFs but no outgoing PDFs are allowed.
 * A more detailed description is available in the dissertation of N. Thuerey, 2007, section 4.4 and table 4.1.
 *
 * neighIt is an iterator pointing to a pre-computed flag field that contains the bitwise OR'ed neighborhood flags of
 * the current cell. See free_surface::getOredNeighborhood for more information.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ConstScalarIt_T, typename ConstPdfIt_T,
          typename ConstFlagIt_T, typename ConstNeighIt_T, typename FlagInfo_T >
real_t advectMass(const FlagField_T* flagField, const ConstScalarIt_T& fillSrc, const ConstPdfIt_T& pdfFieldIt,
                  const ConstFlagIt_T& flagFieldIt, const ConstNeighIt_T& neighIt, const FlagInfo_T& flagInfo,
                  bool useSimpleMassExchange)
{
   using flag_c_t = typename std::remove_const< typename ConstFlagIt_T::value_type >::type;
   using flag_n_t = typename std::remove_const< typename ConstNeighIt_T::value_type >::type;
   using flag_i_t = typename std::remove_const< typename FlagField_T::value_type >::type;

   static_assert(std::is_same< flag_i_t, flag_c_t >::value && std::is_same< flag_i_t, flag_n_t >::value,
                 "Flag types have to be equal.");

   WALBERLA_ASSERT(flagInfo.isInterface(flagFieldIt), "Function advectMass() must only be called for interface cell.");

   // determine the type of the current interface cell (similar to section 4.4 in the dissertation of N. Thuerey,
   // 2007) neighIt is the pre-computed bitwise OR-ed neighborhood of the cell
   bool localNoGasNeig   = !flagInfo.isGas(*neighIt); // this cell has no gas neighbor (should be converted to liquid)
   bool localNoFluidNeig = !flagInfo.isLiquid(*neighIt); // this cell has no fluid neighbor (should be converted to gas)
   bool localStandardCell = !(localNoGasNeig || localNoFluidNeig); // this cell has both gas and fluid neighbors

   // evaluate flag of this cell (flagFieldIt) and not the neighborhood flags (neighIt)
   bool localWettingCell = flagInfo.isKeepInterfaceForWetting(*flagFieldIt); // this cell should be kept for wetting

   if (localNoFluidNeig && localNoGasNeig &&
       !localWettingCell) // this cell has only interface neighbors (interface layer of 3 cells width)
   {
      // WALBERLA_LOG_WARNING("Interface layer of 3 cells width at cell " << fillSrc.cell());
      // set this cell to standard for enabling regular mass exchange
      localNoGasNeig    = false;
      localNoFluidNeig  = false;
      localStandardCell = true;
   }

   real_t deltaMass = real_c(0);
   for (auto dir = LatticeModel_T::Stencil::beginNoCenter(); dir != LatticeModel_T::Stencil::end(); ++dir)
   {
      flag_c_t neighFlag = flagFieldIt.neighbor(*dir);

      bool isFreeSlip = false; // indicates whether dir points to a free slip cell

      // from the viewpoint of the current cell, direction where the PDFs are actually coming from when there is a free
      // slip cell in direction dir; explicitly not a Cell object to emphasize that this denotes a direction
      Vector3< cell_idx_t > freeSlipDir;

      // determine type of cell where the mirrored PDFs at free slip boundary are coming from
      if (flagInfo.isFreeSlip(neighFlag))
      {
         // REMARK: the following implementation is based on lbm/boundary/FreeSlip.h

         // get components of inverse direction of dir
         const int ix = stencil::cx[stencil::inverseDir[*dir]];
         const int iy = stencil::cy[stencil::inverseDir[*dir]];
         const int iz = stencil::cz[stencil::inverseDir[*dir]];

         int wnx = 0; // compute "normal" vector of free slip wall
         int wny = 0;
         int wnz = 0;

         // from the current cell, go into neighboring cell in direction dir and from there determine the type of the
         // neighboring cell in ix, iy and iz direction
         const auto flagFieldFreeSlipPtrX = typename FlagField_T::ConstPtr(
            *flagField, flagFieldIt.x() + dir.cx() + ix, flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
         if (flagInfo.isLiquid(*flagFieldFreeSlipPtrX) || flagInfo.isInterface(*flagFieldFreeSlipPtrX)) { wnx = ix; }

         const auto flagFieldFreeSlipPtrY = typename FlagField_T::ConstPtr(
            *flagField, flagFieldIt.x() + dir.cx(), flagFieldIt.y() + dir.cy() + iy, flagFieldIt.z() + dir.cz());
         if (flagInfo.isLiquid(*flagFieldFreeSlipPtrY) || flagInfo.isInterface(*flagFieldFreeSlipPtrY)) { wny = iy; }

         const auto flagFieldFreeSlipPtrZ = typename FlagField_T::ConstPtr(
            *flagField, flagFieldIt.x() + dir.cx(), flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz() + iz);
         if (flagInfo.isLiquid(*flagFieldFreeSlipPtrZ) || flagInfo.isInterface(*flagFieldFreeSlipPtrZ)) { wnz = iz; }

         // flagFieldFreeSlipPtr denotes the cell from which the PDF is coming from
         const auto flagFieldFreeSlipPtr =
            typename FlagField_T::ConstPtr(*flagField, flagFieldIt.x() + dir.cx() + wnx,
                                           flagFieldIt.y() + dir.cy() + wny, flagFieldIt.z() + dir.cz() + wnz);

         // no mass must be exchanged if:
         // - PDFs are not coming from liquid or interface cells
         // - PDFs are mirrored from this cell (deltaPdf=0)
         if ((!flagInfo.isLiquid(*flagFieldFreeSlipPtr) && !flagInfo.isInterface(*flagFieldFreeSlipPtr)) ||
             flagFieldFreeSlipPtr.cell() == flagFieldIt.cell())
         {
            continue;
         }
         else
         {
            // update neighFlag such that it does contain the flags of the cell that mirrored at the free slip boundary;
            // PDFs from the boundary cell can be used because they were correctly updated by the boundary handling
            neighFlag = *flagFieldFreeSlipPtr;

            // direction of the cell that is mirrored at the free slip boundary
            freeSlipDir = Vector3< cell_idx_t >(cell_idx_c(dir.cx() + wnx), cell_idx_c(dir.cy() + wny),
                                                cell_idx_c(dir.cz() + wnz));
            isFreeSlip  = true;
         }
      }

      // no mass exchange with gas and obstacle cells that are neither inflow nor outflow
      if (flagInfo.isGas(neighFlag) ||
          (flagInfo.isObstacle(neighFlag) && !flagInfo.isInflow(neighFlag) && !flagInfo.isOutflow(neighFlag)))
      {
         continue;
      }

      // PDF pointing from neighbor to current cell
      const real_t neighborPdf = pdfFieldIt.neighbor(*dir, dir.toInvIdx());

      // PDF pointing to neighbor
      const real_t localPdf = pdfFieldIt.getF(dir.toIdx());

      // mass exchange with liquid cells (inflow cells are considered to be liquid)
      if (flagInfo.isLiquid(neighFlag) || flagInfo.isInflow(neighFlag))
      {
         // mass exchange is difference between incoming and outgoing PDFs (see equation (9) in Koerner et al., 2005)
         deltaMass += neighborPdf - localPdf;
         continue;
      }

      // assert cells that are neither gas, obstacle nor interface
      WALBERLA_ASSERT(flagInfo.isInterface(neighFlag) || flagInfo.isOutflow(neighFlag),
                      "In cell " << fillSrc.cell() << ", flag of neighboring cell "
                                 << Cell(fillSrc.x() + dir.cx(), fillSrc.y() + dir.cy(), fillSrc.z() + dir.cz())
                                 << " is not plausible.");

      // direction of the cell from which the PDFs are coming from
      const Vector3< cell_idx_t > relevantDir =
         isFreeSlip ? freeSlipDir :
                      Vector3< cell_idx_t >(cell_idx_c(dir.cx()), cell_idx_c(dir.cy()), cell_idx_c(dir.cz()));

      // determine the type of the neighboring cell (similar to section 4.4 in the dissertation of N. Thuerey, 2007)
      bool neighborNoGasNeig    = !flagInfo.isGas(neighIt.neighbor(relevantDir[0], relevantDir[1], relevantDir[2]));
      bool neighborNoFluidNeig  = !flagInfo.isLiquid(neighIt.neighbor(relevantDir[0], relevantDir[1], relevantDir[2]));
      bool neighborStandardCell = !(neighborNoGasNeig || neighborNoFluidNeig);
      bool neighborWettingCell  = flagInfo.isKeepInterfaceForWetting(flagFieldIt.neighbor(
          relevantDir[0], relevantDir[1],
          relevantDir[2])); // evaluate flag of this cell (flagFieldIt) and not the neighborhood flags (neighIt)

      if (neighborNoGasNeig && neighborNoFluidNeig &&
          !neighborWettingCell) // neighboring cell has only interface neighbors
      {
         // WALBERLA_LOG_WARNING("Interface layer of 3 cells width at cell " << fillSrc.cell());
         //  set neighboring cell to standard for enabling regular mass exchange
         neighborNoGasNeig    = false;
         neighborNoFluidNeig  = false;
         neighborStandardCell = true;
      }

      const real_t localFill    = *fillSrc;
      const real_t neighborFill = fillSrc.neighbor(relevantDir[0], relevantDir[1], relevantDir[2]);

      real_t fillAvg  = real_c(0);
      real_t deltaPdf = real_c(0); // deltaMass = fillAvg * deltaPdf (see equation (9) in Koerner et al., 2005)

      // both cells are interface cells (standard mass exchange)
      // see paper of C. Koerner et al., 2005, equation (9)
      // see dissertation of N. Thuerey, 2007, table 4.1: (standard at x) <-> (standard at x_nb);
      if (useSimpleMassExchange || (localStandardCell && (neighborStandardCell || neighborWettingCell)) ||
          (neighborStandardCell && (localStandardCell || localWettingCell)) || flagInfo.isOutflow(neighFlag))
      {
         if (flagInfo.isOutflow(neighFlag))
         {
            fillAvg = localFill; // use local fill level only, since outflow cells do not have a meaningful fill level
         }
         else { fillAvg = real_c(0.5) * (neighborFill + localFill); }

         deltaPdf = neighborPdf - localPdf;
      }
      else
      {
         // see dissertation of N. Thuerey, 2007, table 4.1:
         //    (standard at x) <-> (no empty neighbors at x_nb)
         //    (no fluid neighbors at x) <-> (standard cell at x_nb)
         //    (no fluid neighbors at x) <-> (no empty neighbors at x_nb)
         // => push local, i.e., this cell empty (if it is not a cell needed for wetting)
         if (((localStandardCell && neighborNoGasNeig) || (localNoFluidNeig && !neighborNoFluidNeig)) &&
             !localWettingCell)
         {
            fillAvg  = real_c(0.5) * (neighborFill + localFill);
            deltaPdf = -localPdf;
         }
         else
         {
            // see dissertation of N. Thuerey, 2007, table 4.1:
            //    (standard at x) <-> (no fluid neighbors at x_nb)
            //    (no empty neighbors at x) <-> (standard cell at x_nb)
            //    (no empty neighbors at x) <-> (no fluid neighbors at x_nb)
            // => push neighboring cell empty (if it is not a cell needed for wetting)
            if (((localStandardCell && neighborNoFluidNeig) || (localNoGasNeig && !neighborNoGasNeig)) &&
                !neighborWettingCell)
            {
               fillAvg  = real_c(0.5) * (neighborFill + localFill);
               deltaPdf = neighborPdf;
            }
            else
            {
               // see dissertation of N. Thuerey, 2007, table 4.1:
               //    (no fluid neighbors at x) <-> (no fluid neighbors at x_nb)
               //    (no empty neighbors at x) <-> (no empty neighbors at x_nb)
               if ((localNoFluidNeig && neighborNoFluidNeig) || (localNoGasNeig && neighborNoGasNeig))
               {
                  fillAvg  = real_c(0.5) * (neighborFill + localFill);
                  deltaPdf = (neighborPdf - localPdf);
               }
               else
               {
                  // treat remaining cases that were not covered above and include wetting cells
                  if (localWettingCell || neighborWettingCell)
                  {
                     fillAvg  = real_c(0.5) * (neighborFill + localFill);
                     deltaPdf = (neighborPdf - localPdf);
                  }
                  else
                  {
                     WALBERLA_ABORT("Unknown mass advection combination of flags (loc=" << *flagFieldIt << ", neig="
                                                                                        << neighFlag << ")");
                  }
               }
            }
         }
      }
      // this cell's deltaMass is the sum over all stencil directions (see dissertation of N. Thuerey, 2007, equation
      // (4.4))
      deltaMass += fillAvg * deltaPdf;
   }

   return deltaMass;
}

} // namespace free_surface
} // namespace walberla
