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
//! \file ForceDensitySweep.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Weight force in interface cells with fill level and density (equation 15 in Koerner et al., 2005).
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Update the force density field as in equation 15 in Koerner et al., 2005. with "acceleration * density * fill level".
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename VectorField_T, typename ScalarField_T >
class ForceDensitySweep
{
 public:
   ForceDensitySweep(BlockDataID forceDensityFieldID, ConstBlockDataID pdfFieldID, ConstBlockDataID flagFieldID,
                     ConstBlockDataID fillFieldID, const FlagInfo< FlagField_T >& flagInfo,
                     const Vector3< real_t >& globalAcceleration)
      : forceDensityFieldID_(forceDensityFieldID), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID),
        fillFieldID_(fillFieldID), flagInfo_(flagInfo), globalAcceleration_(globalAcceleration)
   {}

   void operator()(IBlock* const block)
   {
      using PdfField_T = lbm::PdfField< LatticeModel_T >;

      VectorField_T* const forceDensityField = block->getData< VectorField_T >(forceDensityFieldID_);
      const PdfField_T* const pdfField       = block->getData< const PdfField_T >(pdfFieldID_);
      const FlagField_T* const flagField     = block->getData< const FlagField_T >(flagFieldID_);
      const ScalarField_T* const fillField   = block->getData< const ScalarField_T >(fillFieldID_);

      WALBERLA_FOR_ALL_CELLS(forceDensityFieldIt, forceDensityField, pdfFieldIt, pdfField, flagFieldIt, flagField,
                             fillFieldIt, fillField, {
                                flag_t flag = *flagFieldIt;

                                // set force density in cells to acceleration * density * fillLevel (see equation 15
                                // in Koerner et al., 2005);
                                if (flagInfo_.isInterface(flag))
                                {
                                   const real_t density = pdfField->getDensity(pdfFieldIt.cell());
                                   *forceDensityFieldIt = globalAcceleration_ * *fillFieldIt * density;
                                }

                                else
                                {
                                   if (flagInfo_.isLiquid(flag))
                                   {
                                      const real_t density = pdfField->getDensity(pdfFieldIt.cell());
                                      *forceDensityFieldIt = globalAcceleration_ * density;
                                   }
                                }
                             }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   using flag_t = typename FlagField_T::flag_t;

   BlockDataID forceDensityFieldID_;
   ConstBlockDataID pdfFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID fillFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
   Vector3< real_t > globalAcceleration_;
}; // class ForceDensitySweep

/***********************************************************************************************************************
 * Update the force density field as in equation 15 in Koerner et al., 2005. with "acceleration * density * fill level".
 * Differs from the version above by using a flattened vector field (GhostLayerField< real_t, 3 >). This is necessary
 * because Pystencils does not support VectorField_T (GhostLayerField< Vector3<real_t>, 1 >).
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename VectorFieldFlattened_T, typename ScalarField_T >
class ForceDensityCodegenSweep
{
 public:
   ForceDensityCodegenSweep(BlockDataID forceDensityFieldID, ConstBlockDataID pdfFieldID, ConstBlockDataID flagFieldID,
                            ConstBlockDataID fillFieldID, const FlagInfo< FlagField_T >& flagInfo,
                            const Vector3< real_t >& globalAcceleration)
      : forceDensityFieldID_(forceDensityFieldID), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID),
        fillFieldID_(fillFieldID), flagInfo_(flagInfo), globalAcceleration_(globalAcceleration)
   {}

   void operator()(IBlock* const block)
   {
      using PdfField_T = lbm::PdfField< LatticeModel_T >;

      VectorFieldFlattened_T* const forceDensityField = block->getData< VectorFieldFlattened_T >(forceDensityFieldID_);
      const PdfField_T* const pdfField                = block->getData< const PdfField_T >(pdfFieldID_);
      const FlagField_T* const flagField              = block->getData< const FlagField_T >(flagFieldID_);
      const ScalarField_T* const fillField            = block->getData< const ScalarField_T >(fillFieldID_);

      WALBERLA_FOR_ALL_CELLS(forceDensityFieldIt, forceDensityField, pdfFieldIt, pdfField, flagFieldIt, flagField,
                             fillFieldIt, fillField, {
                                flag_t flag = *flagFieldIt;

                                // set force density in cells to acceleration * density * fillLevel (see equation 15
                                // in Koerner et al., 2005);
                                if (flagInfo_.isInterface(flag))
                                {
                                   const real_t density   = pdfField->getDensity(pdfFieldIt.cell());
                                   forceDensityFieldIt[0] = globalAcceleration_[0] * *fillFieldIt * density;
                                   forceDensityFieldIt[1] = globalAcceleration_[1] * *fillFieldIt * density;
                                   forceDensityFieldIt[2] = globalAcceleration_[2] * *fillFieldIt * density;
                                }
                                else
                                {
                                   if (flagInfo_.isLiquid(flag))
                                   {
                                      const real_t density   = pdfField->getDensity(pdfFieldIt.cell());
                                      forceDensityFieldIt[0] = globalAcceleration_[0] * density;
                                      forceDensityFieldIt[1] = globalAcceleration_[1] * density;
                                      forceDensityFieldIt[2] = globalAcceleration_[2] * density;
                                   }
                                }
                             }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   using flag_t = typename FlagField_T::flag_t;

   BlockDataID forceDensityFieldID_;
   ConstBlockDataID pdfFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID fillFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
   Vector3< real_t > globalAcceleration_;
}; // class ForceDensitySweep

} // namespace free_surface
} // namespace walberla
