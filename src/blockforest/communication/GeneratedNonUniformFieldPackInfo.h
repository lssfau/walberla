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
//! \file GeneratedNonUniformFieldPackInfo.h
//! \ingroup blockforest/communication
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "blockforest/communication/NonUniformPackInfo.h"
#include "field/refinement/PackInfoHelper.h"



namespace walberla {

/*******************************************************************************************************************
 * \brief  Class for packing non-uniform fields, which are not PDF fields, such as velocity or pressure fields.
 *
 * It overrides all needed function for non-uniform packing from the blockforest::communication::NonUniformPackInfo.
 *
 * The functions are gathering the relevant fields and defining the pack and unpack intervals for the packing kernels.
 * The actual packing kernels are generated and imported as template parameter PackingKernels_T.
 */
//*******************************************************************************************************************

template< typename Field_T, typename PackingKernels_T >
class GeneratedNonUniformFieldPackInfo : public blockforest::communication::NonUniformPackInfo
{
 public:
   using value_type = typename Field_T::value_type;

   GeneratedNonUniformFieldPackInfo( const BlockDataID fieldID )
      : fieldID_(fieldID) {};

   ~GeneratedNonUniformFieldPackInfo() = default;

   bool constantDataExchange() const override { return true; };
   bool threadsafeReceiving() const override { return false; };

   void packDataEqualLevelImpl(const Block* sender, stencil::Direction dir, mpi::SendBuffer& buffer) const override;
   void unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer& buffer) override;
   void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir) override;

   void packDataCoarseToFineImpl(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir, mpi::SendBuffer& buffer) const override;
   void unpackDataCoarseToFine      (       Block * fineReceiver, const BlockID & coarseSender, stencil::Direction dir, mpi::RecvBuffer& buffer) override;
   void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir) override;

   void packDataFineToCoarseImpl(const Block* fineSender, const BlockID& coarseReceiver, stencil::Direction dir, mpi::SendBuffer& buffer) const override;
   void unpackDataFineToCoarse      (       Block * coarseReceiver, const BlockID & fineSender,     stencil::Direction dir, mpi::RecvBuffer& buffer) override;
   void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir) override;

   uint_t sizeEqualLevelSend( const Block * sender, stencil::Direction dir) const;
   uint_t sizeCoarseToFineSend ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir) const;
   uint_t sizeCoarseToFineReceive ( Block* fineReceiver, stencil::Direction dir) const;
   uint_t sizeFineToCoarseSend ( const Block * fineSender, stencil::Direction dir) const;

 private:

   const BlockDataID fieldID_;
   PackingKernels_T kernels_;

};
} //namespace walberla

#include "blockforest/communication/GeneratedNonUniformFieldPackInfo.impl.h"
