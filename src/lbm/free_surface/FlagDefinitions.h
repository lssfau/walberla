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
//! \file FlagInfo.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Define free surface flags (e.g. conversion flags).
//
//======================================================================================================================

#include "field/FlagUID.h"

namespace walberla
{
namespace free_surface
{
namespace flagIDs
{
/***********************************************************************************************************************
 * Definition of free surface flag IDs.
 **********************************************************************************************************************/
const field::FlagUID interfaceFlagID                   = field::FlagUID("interface");
const field::FlagUID liquidFlagID                      = field::FlagUID("liquid");
const field::FlagUID gasFlagID                         = field::FlagUID("gas");
const field::FlagUID convertedFlagID                   = field::FlagUID("converted");
const field::FlagUID convertToGasFlagID                = field::FlagUID("convert to gas");
const field::FlagUID convertToLiquidFlagID             = field::FlagUID("convert to liquid");
const field::FlagUID convertedFromGasToInterfaceFlagID = field::FlagUID("convert from gas to interface");
const field::FlagUID keepInterfaceForWettingFlagID     = field::FlagUID("convert to and keep interface for wetting");
const field::FlagUID convertToInterfaceForInflowFlagID = field::FlagUID("convert to interface for inflow");

const Set< field::FlagUID > liquidInterfaceFlagIDs = setUnion< field::FlagUID >(liquidFlagID, interfaceFlagID);

const Set< field::FlagUID > liquidInterfaceGasFlagIDs =
   setUnion(liquidInterfaceFlagIDs, Set< field::FlagUID >(gasFlagID));

} // namespace flagIDs
} // namespace free_surface
} // namespace walberla