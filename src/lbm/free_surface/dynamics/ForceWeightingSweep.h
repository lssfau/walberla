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
//! \file ForceWeightingSweep.h
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
 * Weight the specified force in interface cells according to their density and fill level, as in equation 15 in Koerner
 * et al., 2005.
 * In liquid cells, the force is set to the specified constant global force.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T, typename VectorField_T, typename ScalarField_T >
class ForceWeightingSweep
{
 public:
   ForceWeightingSweep(BlockDataID forceFieldID, ConstBlockDataID pdfFieldID, ConstBlockDataID flagFieldID,
                       ConstBlockDataID fillFieldID, const FlagInfo< FlagField_T >& flagInfo,
                       const Vector3< real_t >& globalForce)
      : forceFieldID_(forceFieldID), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID),
        flagInfo_(flagInfo), globalForce_(globalForce)
   {}

   void operator()(IBlock* const block)
   {
      using PdfField_T = lbm::PdfField< LatticeModel_T >;

      VectorField_T* const forceField      = block->getData< VectorField_T >(forceFieldID_);
      const PdfField_T* const pdfField     = block->getData< const PdfField_T >(pdfFieldID_);
      const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
      const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

      WALBERLA_FOR_ALL_CELLS(forceFieldIt, forceField, pdfFieldIt, pdfField, flagFieldIt, flagField, fillFieldIt,
                             fillField, {
                                flag_t flag = *flagFieldIt;

                                // set force in interface cells to globalForce_ * density * fillLevel (see equation 15
                                // in Koerner et al., 2005)
                                if (flagInfo_.isInterface(flag))
                                {
                                   const real_t density   = pdfField->getDensity(pdfFieldIt.cell());
                                   const real_t fillLevel = *fillFieldIt;
                                   *forceFieldIt          = globalForce_ * fillLevel * density;
                                }
                                else
                                {
                                   // set force to globalForce_ in all non-interface cells
                                   *forceFieldIt = globalForce_;
                                }
                             }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   using flag_t = typename FlagField_T::flag_t;

   BlockDataID forceFieldID_;
   ConstBlockDataID pdfFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID fillFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
   Vector3< real_t > globalForce_;
}; // class ForceWeightingSweep

} // namespace free_surface
} // namespace walberla
