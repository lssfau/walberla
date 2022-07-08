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
//! \file FlagInfo.impl.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Define and manage free surface flags (e.g. conversion flags).
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/mpi/Broadcast.h"

#include "field/FlagField.h"

#include "FlagInfo.h"

namespace walberla
{
namespace free_surface
{
template< typename FlagField_T >
FlagInfo< FlagField_T >::FlagInfo()
   : interfaceFlag(0), liquidFlag(0), gasFlag(0), convertedFlag(0), convertToGasFlag(0), convertToLiquidFlag(0),
     convertFromGasToInterfaceFlag(0), keepInterfaceForWettingFlag(0), convertToInterfaceForInflowFlag(0),
     obstacleFlagMask(0), outflowFlagMask(0), inflowFlagMask(0), freeSlipFlagMask(0)
{}

template< typename FlagField_T >
FlagInfo< FlagField_T >::FlagInfo(const FlagField_T* field, const Set< field::FlagUID >& obstacleIDs,
                                  const Set< field::FlagUID >& outflowIDs, const Set< field::FlagUID >& inflowIDs,
                                  const Set< field::FlagUID >& freeSlipIDs)
   : interfaceFlag(field->getFlag(flagIDs::interfaceFlagID)), liquidFlag(field->getFlag(flagIDs::liquidFlagID)),
     gasFlag(field->getFlag(flagIDs::gasFlagID)), convertedFlag(field->getFlag(flagIDs::convertedFlagID)),
     convertToGasFlag(field->getFlag(flagIDs::convertToGasFlagID)),
     convertToLiquidFlag(field->getFlag(flagIDs::convertToLiquidFlagID)),
     convertFromGasToInterfaceFlag(field->getFlag(flagIDs::convertedFromGasToInterfaceFlagID)),
     keepInterfaceForWettingFlag(field->getFlag(flagIDs::keepInterfaceForWettingFlagID)),
     convertToInterfaceForInflowFlag(field->getFlag(flagIDs::convertToInterfaceForInflowFlagID)),
     obstacleIDs_(obstacleIDs), outflowIDs_(outflowIDs), inflowIDs_(inflowIDs), freeSlipIDs_(freeSlipIDs)
{
   // create obstacleFlagMask from obstacleIDs using bitwise OR
   obstacleFlagMask = 0;
   for (auto obstacleID = obstacleIDs.begin(); obstacleID != obstacleIDs.end(); ++obstacleID)
   {
      obstacleFlagMask = flag_t(obstacleFlagMask | field->getFlag(*obstacleID));
   }

   // create outflowFlagMask from outflowIDs using bitwise OR
   outflowFlagMask = 0;
   for (auto outflowID = outflowIDs.begin(); outflowID != outflowIDs.end(); ++outflowID)
   {
      outflowFlagMask = flag_t(outflowFlagMask | field->getFlag(*outflowID));
   }

   // create inflowFlagMask from inflowIDs using bitwise OR
   inflowFlagMask = 0;
   for (auto inflowID = inflowIDs.begin(); inflowID != inflowIDs.end(); ++inflowID)
   {
      inflowFlagMask = flag_t(inflowFlagMask | field->getFlag(*inflowID));
   }

   // create freeSlipFlagMask from freeSlipIDs using bitwise OR
   freeSlipFlagMask = 0;
   for (auto freeSlipID = freeSlipIDs.begin(); freeSlipID != freeSlipIDs.end(); ++freeSlipID)
   {
      freeSlipFlagMask = flag_t(freeSlipFlagMask | field->getFlag(*freeSlipID));
   }
}

template< typename FlagField_T >
FlagInfo< FlagField_T >::FlagInfo(const Set< field::FlagUID >& obstacleIDs, const Set< field::FlagUID >& outflowIDs,
                                  const Set< field::FlagUID >& inflowIDs, const Set< field::FlagUID >& freeSlipIDs)
   : obstacleIDs_(obstacleIDs), outflowIDs_(outflowIDs), inflowIDs_(inflowIDs), freeSlipIDs_(freeSlipIDs)
{
   // define flags
   flag_t nextFreeBit = flag_t(0);
   interfaceFlag = flag_t(flag_t(1) << nextFreeBit++); // outermost flag_t is necessary to avoid warning C4334 on MSVC
   liquidFlag    = flag_t(flag_t(1) << nextFreeBit++);
   gasFlag       = flag_t(flag_t(1) << nextFreeBit++);
   convertedFlag = flag_t(flag_t(1) << nextFreeBit++);
   convertToGasFlag                = flag_t(flag_t(1) << nextFreeBit++);
   convertToLiquidFlag             = flag_t(flag_t(1) << nextFreeBit++);
   convertFromGasToInterfaceFlag   = flag_t(flag_t(1) << nextFreeBit++);
   keepInterfaceForWettingFlag     = flag_t(flag_t(1) << nextFreeBit++);
   convertToInterfaceForInflowFlag = flag_t(flag_t(1) << nextFreeBit++);

   obstacleFlagMask = flag_t(0);
   outflowFlagMask  = flag_t(0);
   inflowFlagMask   = flag_t(0);
   freeSlipFlagMask = flag_t(0);

   // define flags for obstacles, outflow, inflow and freeSlip
   auto setUnion = obstacleIDs + outflowIDs + inflowIDs + freeSlipIDs;
   for (auto id = setUnion.begin(); id != setUnion.end(); ++id)
   {
      if (obstacleIDs.contains(*id)) { obstacleFlagMask = obstacleFlagMask | flag_t((flag_t(1) << nextFreeBit)); }

      if (outflowIDs.contains(*id)) { outflowFlagMask = outflowFlagMask | flag_t((flag_t(1) << nextFreeBit)); }

      if (inflowIDs.contains(*id)) { inflowFlagMask = inflowFlagMask | flag_t((flag_t(1) << nextFreeBit)); }

      if (freeSlipIDs.contains(*id)) { freeSlipFlagMask = freeSlipFlagMask | flag_t((flag_t(1) << nextFreeBit)); }

      nextFreeBit++;
   }
}

template< typename FlagField_T >
void FlagInfo< FlagField_T >::registerFlags(FlagField_T* field, const Set< field::FlagUID >& obstacleIDs,
                                            const Set< field::FlagUID >& outflowIDs,
                                            const Set< field::FlagUID >& inflowIDs,
                                            const Set< field::FlagUID >& freeSlipIDs)
{
   flag_t nextFreeBit = flag_t(0);
   field->registerFlag(flagIDs::interfaceFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::liquidFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::gasFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::convertedFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::convertToGasFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::convertToLiquidFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::convertedFromGasToInterfaceFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::keepInterfaceForWettingFlagID, nextFreeBit++);
   field->registerFlag(flagIDs::convertToInterfaceForInflowFlagID, nextFreeBit++);

   // extract flags
   auto setUnion = obstacleIDs + outflowIDs + inflowIDs + freeSlipIDs;
   for (auto id = setUnion.begin(); id != setUnion.end(); ++id)
   {
      field->registerFlag(*id, nextFreeBit++);
   }
}

template< typename FlagField_T >
void FlagInfo< FlagField_T >::registerFlags(FlagField_T* field) const
{
   registerFlags(field, obstacleIDs_, outflowIDs_, inflowIDs_, freeSlipIDs_);
}

template< typename FlagField_T >
bool FlagInfo< FlagField_T >::isConsistentAcrossBlocksAndProcesses(
   const std::weak_ptr< StructuredBlockStorage >& blockStoragePtr, ConstBlockDataID flagFieldID) const
{
   // check consistency across processes
   FlagInfo rootFlagInfo = *this;

   // root broadcasts its FlagInfo to all other processes
   mpi::broadcastObject(rootFlagInfo);

   // this process' FlagInfo is not identical to the one of root
   if (rootFlagInfo != *this) { return false; }

   auto blockStorage = blockStoragePtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockStorage);

   // check consistency across blocks
   for (auto blockIt = blockStorage->begin(); blockIt != blockStorage->end(); ++blockIt)
   {
      const FlagField_T* const flagField = blockIt->getData< const FlagField_T >(flagFieldID);
      FlagInfo fi(flagField, obstacleIDs_, outflowIDs_, inflowIDs_, freeSlipIDs_);
      if (fi != *this) { return false; }
   }

   return true;
}

template< typename FlagField_T >
bool FlagInfo< FlagField_T >::operator==(const FlagInfo& o) const
{
   return interfaceFlag == o.interfaceFlag && gasFlag == o.gasFlag && liquidFlag == o.liquidFlag &&
          convertToGasFlag == o.convertToGasFlag && convertToLiquidFlag == o.convertToLiquidFlag &&
          convertFromGasToInterfaceFlag == o.convertFromGasToInterfaceFlag &&
          keepInterfaceForWettingFlag == o.keepInterfaceForWettingFlag &&
          convertToInterfaceForInflowFlag == o.convertToInterfaceForInflowFlag &&
          obstacleFlagMask == o.obstacleFlagMask && outflowFlagMask == o.outflowFlagMask &&
          inflowFlagMask == o.inflowFlagMask && freeSlipFlagMask == o.freeSlipFlagMask;
}

template< typename FlagField_T >
bool FlagInfo< FlagField_T >::operator!=(const FlagInfo& o) const
{
   return !(*this == o);
}

} // namespace free_surface
} // namespace walberla