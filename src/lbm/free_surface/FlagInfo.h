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
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Manage free surface flags (e.g. conversion flags).
//
//======================================================================================================================

#pragma once

#include "core/mpi/Broadcast.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

#include "FlagDefinitions.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Manage free surface flags (e.g. conversion flags).
 **********************************************************************************************************************/
template< typename FlagField_T >
class FlagInfo
{
 public:
   using flag_t = typename FlagField_T::flag_t;

   FlagInfo();
   FlagInfo(const Set< field::FlagUID >& obstacleIDs,
            const Set< field::FlagUID >& outflowIDs  = Set< field::FlagUID >::emptySet(),
            const Set< field::FlagUID >& inflowIDs   = Set< field::FlagUID >::emptySet(),
            const Set< field::FlagUID >& freeSlipIDs = Set< field::FlagUID >::emptySet());

   bool operator==(const FlagInfo& o) const;
   bool operator!=(const FlagInfo& o) const;

   flag_t interfaceFlag;
   flag_t liquidFlag;
   flag_t gasFlag;

   flag_t convertedFlag;                   // interface cell that is already converted
   flag_t convertToGasFlag;                // interface cell with too low fill level, should be converted
   flag_t convertToLiquidFlag;             // interface cell with too high fill level, should be converted
   flag_t convertFromGasToInterfaceFlag;   // interface cell that was a gas cell and needs refilling of its pdfs
   flag_t keepInterfaceForWettingFlag;     // gas/liquid cell that needs to be converted to interface cell for a smooth
                                           // continuation of the wetting surface (see dissertation of S. Donath, 2011,
                                           // section 6.3.5.3)
   flag_t convertToInterfaceForInflowFlag; // gas cell that needs to be converted to interface cell to enable inflow

   flag_t obstacleFlagMask;
   flag_t outflowFlagMask;
   flag_t inflowFlagMask;
   flag_t freeSlipFlagMask; // free slip obstacle cell (needs to be treated separately since PDFs going from gas into
                            // boundary are not available and must be reconstructed)

   template< typename FieldItOrPtr_T >
   inline bool isInterface(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, interfaceFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool isLiquid(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, liquidFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool isGas(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, gasFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool isObstacle(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isPartOfMaskSet(flagItOrPtr, obstacleFlagMask);
   }
   template< typename FieldItOrPtr_T >
   inline bool isOutflow(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isPartOfMaskSet(flagItOrPtr, outflowFlagMask);
   }
   template< typename FieldItOrPtr_T >
   inline bool isInflow(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isPartOfMaskSet(flagItOrPtr, inflowFlagMask);
   }
   template< typename FieldItOrPtr_T >
   inline bool isFreeSlip(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isPartOfMaskSet(flagItOrPtr, freeSlipFlagMask);
   }
   template< typename FieldItOrPtr_T >
   inline bool hasConverted(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, convertedFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool hasConvertedToGas(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, convertToGasFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool hasConvertedToLiquid(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, convertToLiquidFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool hasConvertedFromGasToInterface(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, convertFromGasToInterfaceFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool isKeepInterfaceForWetting(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, keepInterfaceForWettingFlag);
   }
   template< typename FieldItOrPtr_T >
   inline bool isConvertToInterfaceForInflow(const FieldItOrPtr_T& flagItOrPtr) const
   {
      return isFlagSet(flagItOrPtr, convertToInterfaceForInflowFlag);
   }

   inline bool isInterface(const flag_t val) const { return isFlagSet(val, interfaceFlag); }
   inline bool isLiquid(const flag_t val) const { return isFlagSet(val, liquidFlag); }
   inline bool isGas(const flag_t val) const { return isFlagSet(val, gasFlag); }
   inline bool isObstacle(const flag_t val) const { return isPartOfMaskSet(val, obstacleFlagMask); }
   inline bool isOutflow(const flag_t val) const { return isPartOfMaskSet(val, outflowFlagMask); }
   inline bool isInflow(const flag_t val) const { return isPartOfMaskSet(val, inflowFlagMask); }
   inline bool isFreeSlip(const flag_t val) const { return isPartOfMaskSet(val, freeSlipFlagMask); }
   inline bool isKeepInterfaceForWetting(const flag_t val) const { return isFlagSet(val, keepInterfaceForWettingFlag); }
   inline bool hasConvertedToGas(const flag_t val) const { return isFlagSet(val, convertToGasFlag); }
   inline bool hasConvertedToLiquid(const flag_t val) const { return isFlagSet(val, convertToLiquidFlag); }
   inline bool hasConvertedFromGasToInterface(const flag_t val) const
   {
      return isFlagSet(val, convertFromGasToInterfaceFlag);
   }
   inline bool isConvertToInterfaceForInflow(const flag_t val) const
   {
      return isFlagSet(val, convertToInterfaceForInflowFlag);
   }

   // check whether FlagInfo is identical on all blocks and all processes
   bool isConsistentAcrossBlocksAndProcesses(const std::weak_ptr< StructuredBlockStorage >& blockStorage,
                                             ConstBlockDataID flagFieldID) const;

   // register flags in flag field
   static void registerFlags(FlagField_T* field, const Set< field::FlagUID >& obstacleIDs,
                             const Set< field::FlagUID >& outflowIDs  = Set< field::FlagUID >::emptySet(),
                             const Set< field::FlagUID >& inflowIDs   = Set< field::FlagUID >::emptySet(),
                             const Set< field::FlagUID >& freeSlipIDs = Set< field::FlagUID >::emptySet());

   void registerFlags(FlagField_T* field) const;

   Set< field::FlagUID > getObstacleIDSet() const { return obstacleIDs_; }
   Set< field::FlagUID > getOutflowIDs() const { return outflowIDs_; }
   Set< field::FlagUID > getInflowIDs() const { return inflowIDs_; }
   Set< field::FlagUID > getFreeSlipIDs() const { return freeSlipIDs_; }

 protected:
   FlagInfo(const FlagField_T* field, const Set< field::FlagUID >& obstacleIDs, const Set< field::FlagUID >& outflowIDs,
            const Set< field::FlagUID >& inflowIDs, const Set< field::FlagUID >& freeSlipIDs);

   // create sets of flag IDs with flags from free_surface/boundary/FreeSurfaceBoundaryHandling.impl.h
   Set< field::FlagUID > obstacleIDs_;
   Set< field::FlagUID > outflowIDs_;
   Set< field::FlagUID > inflowIDs_;
   Set< field::FlagUID > freeSlipIDs_;
};

template< typename FlagField_T >
mpi::SendBuffer& operator<<(mpi::SendBuffer& buf, const FlagInfo< FlagField_T >& flagInfo)
{
   buf << flagInfo.interfaceFlag << flagInfo.liquidFlag << flagInfo.gasFlag << flagInfo.convertedFlag
       << flagInfo.convertToGasFlag << flagInfo.convertToLiquidFlag << flagInfo.convertFromGasToInterfaceFlag
       << flagInfo.keepInterfaceForWettingFlag << flagInfo.keepInterfaceForWettingFlag
       << flagInfo.convertToInterfaceForInflowFlag << flagInfo.obstacleFlagMask << flagInfo.outflowFlagMask
       << flagInfo.inflowFlagMask << flagInfo.freeSlipFlagMask;

   return buf;
}

template< typename FlagField_T >
mpi::RecvBuffer& operator>>(mpi::RecvBuffer& buf, FlagInfo< FlagField_T >& flagInfo)
{
   buf >> flagInfo.interfaceFlag >> flagInfo.liquidFlag >> flagInfo.gasFlag >> flagInfo.convertedFlag >>
      flagInfo.convertToGasFlag >> flagInfo.convertToLiquidFlag >> flagInfo.convertFromGasToInterfaceFlag >>
      flagInfo.keepInterfaceForWettingFlag >> flagInfo.keepInterfaceForWettingFlag >>
      flagInfo.convertToInterfaceForInflowFlag >> flagInfo.obstacleFlagMask >> flagInfo.outflowFlagMask >>
      flagInfo.inflowFlagMask >> flagInfo.freeSlipFlagMask;

   return buf;
}

} // namespace free_surface
} // namespace walberla

#include "FlagInfo.impl.h"