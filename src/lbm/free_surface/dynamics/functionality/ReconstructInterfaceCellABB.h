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
//! \file ReconstructInterfaceCellABB.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Free surface boundary condition as in Koerner et al., 2005. Similar to anti-bounce-back pressure condition.
//
//======================================================================================================================

#pragma once

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/dynamics/PdfReconstructionModel.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"

#include "stencil/Directions.h"

#include "GetLaplacePressure.h"

namespace walberla
{
namespace free_surface
{
// get index of largest entry in n_dot_ci with isInterfaceOrLiquid==true && isPdfAvailable==false
uint_t getIndexOfMaximum(const std::vector< bool >& isInterfaceOrLiquid, const std::vector< bool >& isPdfAvailable,
                         const std::vector< real_t >& n_dot_ci);

// get index of smallest entry in n_dot_ci with isInterfaceOrLiquid==true && isPdfAvailable==false
uint_t getIndexOfMinimum(const std::vector< bool >& isInterfaceOrLiquid, const std::vector< bool >& isPdfAvailable,
                         const std::vector< real_t >& n_dot_ci);

// reconstruct PDFs according to pressure anti bounce back boundary condition (page 31, equation 4.5 in dissertation of
// N. Thuerey, 2007)
template< typename LatticeModel_T, typename ConstPdfIt_T >
inline real_t reconstructPressureAntiBounceBack(const stencil::Iterator< typename LatticeModel_T::Stencil >& dir,
                                                const ConstPdfIt_T& pdfFieldIt, const Vector3< real_t >& u,
                                                real_t rhoGas, real_t dir_independent)
{
   const real_t vel = real_c(dir.cx()) * u[0] + real_c(dir.cy()) * u[1] + real_c(dir.cz()) * u[2];

   // compute f^{eq}_i + f^{eq}_{\overline{i}} using rhoGas (without linear terms as they cancel out)
   const real_t tmp =
      real_c(2.0) * LatticeModel_T::w[dir.toIdx()] * rhoGas * (dir_independent + real_c(4.5) * vel * vel);

   // reconstruct PDFs (page 31, equation 4.5 in dissertation of N. Thuerey, 2007)
   return tmp - pdfFieldIt.getF(dir.toIdx());
}

/***********************************************************************************************************************
 * Free surface boundary condition as described in the publication of Koerner et al., 2005. Missing PDFs are
 * reconstructed according to an anti-bounce-back pressure boundary condition at the free surface.
 *
 * Be aware that in Koerner et al., 2005, the PDFs are reconstructed in gas cells (neighboring to interface cells).
 * These PDFs are then streamed into the interface cells in the LBM stream. Here, we directly reconstruct the PDFs in
 * interface cells such that our implementation follows the notation in the dissertation of N. Thuerey, 2007, page 31.
 * ********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename ConstPdfIt_T, typename ConstFlagIt_T,
          typename ConstVectorIt_T, typename OutputArray_T >
void reconstructInterfaceCellLegacy(const FlagField_T* flagField, const ConstPdfIt_T& pdfFieldIt,
                                    const ConstFlagIt_T& flagFieldIt, const ConstVectorIt_T& normalFieldIt,
                                    const FlagInfo< FlagField_T >& flagInfo, const real_t rhoGas, OutputArray_T& f,
                                    const PdfReconstructionModel& pdfReconstructionModel)
{
   using Stencil_T = typename LatticeModel_T::Stencil;

   // get velocity and density in interface cell
   Vector3< real_t > u;
   auto pdfField = dynamic_cast< const lbm::PdfField< LatticeModel_T >* >(pdfFieldIt.getField());
   WALBERLA_ASSERT_NOT_NULLPTR(pdfField);
   const real_t rho = lbm::getDensityAndMomentumDensity(u, pdfField->latticeModel(), pdfFieldIt);
   u /= rho;

   const real_t dir_independent =
      real_c(1.0) - real_c(1.5) * u.sqrLength(); // direction independent value used for PDF reconstruction

   // get type of the model that determines the PDF reconstruction
   const PdfReconstructionModel::ReconstructionModel reconstructionModel = pdfReconstructionModel.getModelType();

   // vector that stores the dot product between interface normal and lattice direction for each lattice direction
   std::vector< real_t > n_dot_ci;

   // vector that stores which index from loop "for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end();
   // ++dir)" is currently available, i.e.:
   // - reconstructed (or scheduled for reconstruction)
   // - coming from boundary condition
   // - available due to fluid or interface neighbor
   std::vector< bool > isPdfAvailable;

   // vector that stores which index from loop "for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end();
   // ++dir)" points to a neighboring interface or fluid cell
   std::vector< bool > isInterfaceOrLiquid;

   // count number of reconstructed links
   uint_t numReconstructed = uint_c(0);

   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      const auto neighborFlag = flagFieldIt.neighbor(*dir);

      if (flagInfo.isObstacle(neighborFlag))
      {
         // free slip boundaries need special treatment because PDFs traveling from gas cells into the free slip
         // boundary must be reconstructed, for instance:
         // [I][G]  with I: interface cell; G: gas cell; f: free slip cell
         // [f][f]
         // During streaming, the interface cell's PDF with direction (-1, 1) is coming from the right free slip cell.
         // For a free slip boundary, this PDF is identical to the PDF with direction (-1, -1) in the gas cell. Since
         // gas-side PDFs are not available, such PDFs must be reconstructed.
         // Non-gas cells do not need to be treated here as they are treated correctly by the boundary handling.

         if (flagInfo.isFreeSlip(neighborFlag))
         {
            // REMARK: the following implementation is based on lbm/boundary/FreeSlip.h

            // get components of inverse direction of dir
            const int ix = stencil::cx[stencil::inverseDir[*dir]];
            const int iy = stencil::cy[stencil::inverseDir[*dir]];
            const int iz = stencil::cz[stencil::inverseDir[*dir]];

            int wnx = 0; // compute "normal" vector of free slip wall
            int wny = 0;
            int wnz = 0;

            // from the current cell, go into neighboring cell in direction dir and from there check whether the
            // neighbors in ix, iy and iz are gas cells
            const auto flagFieldFreeSlipPtrX = typename FlagField_T::ConstPtr(
               *flagField, flagFieldIt.x() + dir.cx() + ix, flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz());
            if (flagInfo.isGas(*flagFieldFreeSlipPtrX)) { wnx = ix; }

            const auto flagFieldFreeSlipPtrY = typename FlagField_T::ConstPtr(
               *flagField, flagFieldIt.x() + dir.cx(), flagFieldIt.y() + dir.cy() + iy, flagFieldIt.z() + dir.cz());
            if (flagInfo.isGas(*flagFieldFreeSlipPtrY)) { wny = iy; }

            const auto flagFieldFreeSlipPtrZ = typename FlagField_T::ConstPtr(
               *flagField, flagFieldIt.x() + dir.cx(), flagFieldIt.y() + dir.cy(), flagFieldIt.z() + dir.cz() + iz);
            if (flagInfo.isGas(*flagFieldFreeSlipPtrZ)) { wnz = iz; }

            if (wnx != 0 || wny != 0 || wnz != 0)
            {
               // boundaryNeighbor denotes the cell from which the PDF is coming from
               const auto flagFieldFreeSlipPtr =
                  typename FlagField_T::ConstPtr(*flagField, flagFieldIt.x() + dir.cx() + wnx,
                                                 flagFieldIt.y() + dir.cy() + wny, flagFieldIt.z() + dir.cz() + wnz);

               if (flagInfo.isGas(*flagFieldFreeSlipPtr))
               {
                  // reconstruct PDF
                  f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(
                     dir, pdfFieldIt, u, rhoGas, dir_independent);
                  isPdfAvailable.push_back(true);
                  isInterfaceOrLiquid.push_back(false);
                  n_dot_ci.push_back(
                     real_c(0)); // dummy entry for having index as in vectors isPDFAvailable and isInterfaceOrLiquid
                  ++numReconstructed;
                  continue;
               }
               // else: do nothing here, i.e., make usual obstacle boundary treatment below
            } // else: concave corner, all surrounding PDFs are known, i.e., make usual obstacle boundary treatment
              // below
         }

         f[dir.toInvIdx()] = pdfFieldIt.neighbor(*dir, dir.toInvIdx()); // use PDFs defined by boundary handling
         isPdfAvailable.push_back(true);
         isInterfaceOrLiquid.push_back(false);
         n_dot_ci.push_back(
            real_c(0)); // dummy entry for having same indices as in vectors isPDFAvailable and isInterfaceOrLiquid
         continue;
      }
      else
      {
         if (flagInfo.isGas(neighborFlag))
         {
            f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(
               dir, pdfFieldIt, u, rhoGas, dir_independent);
            isPdfAvailable.push_back(true);
            isInterfaceOrLiquid.push_back(false);
            n_dot_ci.push_back(
               real_c(0)); // dummy entry for having index as in vectors isPDFAvailable and isInterfaceOrLiquid
            ++numReconstructed;
            continue;
         }
      }

      // dot product between interface normal and lattice direction
      real_t dotProduct = (*normalFieldIt)[0] * real_c(dir.cx()) + (*normalFieldIt)[1] * real_c(dir.cy()) +
                          (*normalFieldIt)[2] * real_c(dir.cz());

      // avoid numerical inaccuracies in the computation of the scalar product n_dot_ci
      if (realIsEqual(dotProduct, real_c(0), real_c(1e-14))) { dotProduct = real_c(0); }

      // approach from Koerner; reconstruct all PDFs in direction opposite to interface normal and center PDF (not
      // stated in the paper explicitly but follows from n*e_i>=0 for i=0)
      if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::NormalBasedReconstructCenter)
      {
         if (dotProduct >= real_c(0))
         {
            f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(
               dir, pdfFieldIt, u, rhoGas, dir_independent);
         }
         else
         {
            // regular LBM stream with PDFs from neighbor
            f[dir.toInvIdx()] = pdfFieldIt.neighbor(*dir, dir.toInvIdx());
         }
         continue;
      }

      if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::NormalBasedKeepCenter)
      {
         if (*dir == stencil::C)
         {
            // use old center PDF
            f[Stencil_T::idx[stencil::C]] = pdfFieldIt[Stencil_T::idx[stencil::C]];
         }
         else
         {
            if (dotProduct >= real_c(0))
            {
               f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(
                  dir, pdfFieldIt, u, rhoGas, dir_independent);
            }
            else
            {
               // regular LBM stream with PDFs from neighbor
               f[dir.toInvIdx()] = pdfFieldIt.neighbor(*dir, dir.toInvIdx());
            }
         }
         continue;
      }

      // reconstruct all non-obstacle PDFs, including those that come from liquid and are already known
      if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::All)
      {
         f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(dir, pdfFieldIt, u,
                                                                                               rhoGas, dir_independent);
         continue;
      }

      // reconstruct only those gas-side PDFs that are really missing
      if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissing)
      {
         // regular LBM stream with PDFs from neighboring interface or liquid cell
         f[dir.toInvIdx()] = pdfFieldIt.neighbor(*dir, dir.toInvIdx());
         continue;
      }

      // reconstruct only those gas-side PDFs that are really missing but make sure that at least a
      // specified number of PDFs are reconstructed (even if available PDFs are overwritten)
      if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin)
      {
         isPdfAvailable.push_back(false); // PDF has not yet been marked for reconstruction
         isInterfaceOrLiquid.push_back(true);
         n_dot_ci.push_back(dotProduct);
         continue;
      }

      WALBERLA_ABORT("Unknown pdfReconstructionModel.")
   }

   if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin)
   {
      WALBERLA_ASSERT_EQUAL(Stencil_T::Size, uint_c(n_dot_ci.size()));
      WALBERLA_ASSERT_EQUAL(Stencil_T::Size, uint_c(isInterfaceOrLiquid.size()));
      WALBERLA_ASSERT_EQUAL(Stencil_T::Size, uint_c(isPdfAvailable.size()));

      const uint_t numMinReconstruct = pdfReconstructionModel.getNumMinReconstruct();

      // number of remaining PDFs that need to be reconstructed according to the specified model (number can be negative
      // => do not use uint_t)
      int numRemainingReconstruct = int_c(numMinReconstruct) - int_c(numReconstructed);

      // count the number of neighboring cells that are interface or liquid (and not obstacle or gas)
      const uint_t numLiquidNeighbors =
         uint_c(std::count_if(isInterfaceOrLiquid.begin(), isInterfaceOrLiquid.end(), [](bool a) { return a; }));

      // // REMARK: this was commented because it regularly occurred in practical simulations (e.g. BreakingDam)
      //    if (numRemainingReconstruct > int_c(0) && numRemainingReconstruct < int_c(numLiquidNeighbors))
      //    {
      //       // There are less neighboring liquid and interface cells than needed to reconstruct PDFs in the
      //       // free surface boundary condition. You have probably specified a large minimum number of PDFs to be
      //       // reconstructed. This happens e.g. near solid boundaries (especially corners). There, the number of
      //       // surrounding non-obstacle cells might be less than the number of PDFs that you want to have
      //       // reconstructed.
      //       WALBERLA_LOG_WARNING_ON_ROOT("Less PDFs reconstructed in cell "
      //                                     << pdfFieldIt.cell()
      //                                     << " than specified in the PDF reconstruction model. See comment in "
      //                                        "source code of ReconstructInterfaceCellABB.h for further information. "
      //                                        "Here, as many PDFs as possible are reconstructed now.");
      //    }

      // count additionally reconstructed PDFs (that come from interface or liquid)
      uint_t numAdditionalReconstructed = uint_c(0);

      // define which PDFs to additionally reconstruct (overwrite known information) according to the specified model
      while (numRemainingReconstruct > int_c(0) && numAdditionalReconstructed < numLiquidNeighbors)
      {
         if (pdfReconstructionModel.getFallbackModel() == PdfReconstructionModel::FallbackModel::Largest)
         {
            // get index of largest n_dot_ci with isInterfaceOrLiquid==true && isPdfAvailable==false
            const uint_t maxIndex    = getIndexOfMaximum(isInterfaceOrLiquid, isPdfAvailable, n_dot_ci);
            isPdfAvailable[maxIndex] = true; // reconstruct this PDF later
            ++numReconstructed;
            ++numAdditionalReconstructed;
         }
         else
         {
            if (pdfReconstructionModel.getFallbackModel() == PdfReconstructionModel::FallbackModel::Smallest)
            {
               // get index of smallest n_dot_ci with isInterfaceOrLiquid==true && isPdfAvailable==false
               const uint_t minIndex    = getIndexOfMinimum(isInterfaceOrLiquid, isPdfAvailable, n_dot_ci);
               isPdfAvailable[minIndex] = true; // reconstruct this PDF later
               ++numReconstructed;
               ++numAdditionalReconstructed;
            }
            else
            {
               // use approach from Koerner
               if (pdfReconstructionModel.getFallbackModel() ==
                   PdfReconstructionModel::FallbackModel::NormalBasedKeepCenter)
               {
                  uint_t index = uint_c(0);
                  for (const real_t& value : n_dot_ci)
                  {
                     if (value >= real_c(0) && index > uint_c(1)) // skip center PDF with index=0
                     {
                        isPdfAvailable[index] = true; // reconstruct this PDF later
                     }

                     ++index;
                  }
                  break; // exit while loop
               }
               else { WALBERLA_ABORT("Unknown fallbackModel in pdfReconstructionModel.") }
            }
         }

         numRemainingReconstruct = int_c(numMinReconstruct) - int_c(numReconstructed);
      }

      // reconstruct additional PDFs
      uint_t index = uint_c(0);
      for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir, ++index)
      {
         const auto neighborFlag = flagFieldIt.neighbor(*dir);

         // skip links pointing to obstacle and gas neighbors; they were treated above already
         if (flagInfo.isObstacle(neighborFlag) || flagInfo.isGas(neighborFlag)) { continue; }
         else
         {
            // reconstruct links that were marked for reconstruction
            if (isPdfAvailable[index] && isInterfaceOrLiquid[index])
            {
               f[dir.toInvIdx()] = reconstructPressureAntiBounceBack< LatticeModel_T, ConstPdfIt_T >(
                  dir, pdfFieldIt, u, rhoGas, dir_independent);
            }
            else
            {
               if (!isPdfAvailable[index] && isInterfaceOrLiquid[index])
               {
                  // regular LBM stream with PDFs from neighbor
                  f[dir.toInvIdx()] = pdfFieldIt.neighbor(*dir, dir.toInvIdx());
                  continue;
               }
               WALBERLA_ABORT("Error in PDF reconstruction. This point should never be reached.")
            }
         }
      }
   }
}

uint_t getIndexOfMaximum(const std::vector< bool >& isInterfaceOrLiquid, const std::vector< bool >& isPdfAvailable,
                         const std::vector< real_t >& n_dot_ci)
{
   real_t maximum = -std::numeric_limits< real_t >::max();
   uint_t index   = std::numeric_limits< uint_t >::max();

   for (uint_t i = uint_c(0); i != isInterfaceOrLiquid.size(); ++i)
   {
      if (isInterfaceOrLiquid[i] && !isPdfAvailable[i])
      {
         const real_t absValue = std::abs(n_dot_ci[i]);
         if (absValue > maximum)
         {
            maximum = absValue;
            index   = i;
         }
      }
   }

   // less Pdfs available for being reconstructed than specified by the user; these assertions should never fail, as the
   // conditionals in reconstructInterfaceCellLegacy() should avoid calling this function
   WALBERLA_ASSERT(maximum > -real_c(std::numeric_limits< real_t >::min()));
   WALBERLA_ASSERT(index != std::numeric_limits< uint_t >::max());

   return index;
}

uint_t getIndexOfMinimum(const std::vector< bool >& isInterfaceOrLiquid, const std::vector< bool >& isPdfAvailable,
                         const std::vector< real_t >& n_dot_ci)
{
   real_t minimum = std::numeric_limits< real_t >::max();
   uint_t index   = std::numeric_limits< uint_t >::max();

   for (uint_t i = uint_c(0); i != isInterfaceOrLiquid.size(); ++i)
   {
      if (isInterfaceOrLiquid[i] && !isPdfAvailable[i])
      {
         const real_t absValue = std::abs(n_dot_ci[i]);
         if (absValue < minimum)
         {
            minimum = absValue;
            index   = i;
         }
      }
   }

   // fewer PDFs available for being reconstructed than specified by the user; these assertions should never fail, as
   // the conditionals in reconstructInterfaceCellLegacy() should avoid calling this function
   WALBERLA_ASSERT(minimum < real_c(std::numeric_limits< real_t >::max()));
   WALBERLA_ASSERT(index != std::numeric_limits< uint_t >::max());

   return index;
}

} // namespace free_surface
} // namespace walberla
