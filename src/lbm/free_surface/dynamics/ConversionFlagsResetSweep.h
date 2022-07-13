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
//! \file ResetFlagSweep.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Reset all free surface flags that mark cell conversions.
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

#include "lbm/free_surface/FlagInfo.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Reset all free surface flags that signal cell conversions. The flag "keepInterfaceForWettingFlag" is explicitly not
 * reset since this flag must persist in the next time step.
 **********************************************************************************************************************/
template< typename FlagField_T >
class ConversionFlagsResetSweep
{
 public:
   ConversionFlagsResetSweep(BlockDataID flagFieldID, const FlagInfo< FlagField_T >& flagInfo)
      : flagFieldID_(flagFieldID), flagInfo_(flagInfo)
   {}

   void operator()(IBlock* const block)
   {
      FlagField_T* const flagField = block->getData< FlagField_T >(flagFieldID_);

      // reset all conversion flags (except flagInfo_.keepInterfaceForWettingFlag)
      const flag_t allConversionFlags = flagInfo_.convertToGasFlag | flagInfo_.convertToLiquidFlag |
                                        flagInfo_.convertedFlag | flagInfo_.convertFromGasToInterfaceFlag |
                                        flagInfo_.convertToInterfaceForInflowFlag;
      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField,
                             { removeMask(flagFieldIt, allConversionFlags); }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   using flag_t = typename FlagField_T::flag_t;

   BlockDataID flagFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
}; // class ConversionFlagsResetSweep

} // namespace free_surface
} // namespace walberla
