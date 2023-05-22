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
//! \file CombinedInPlacePackInfo.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once
#include "communication/UniformPackInfo.h"

namespace walberla::lbm_generated {

template< typename LatticeStorageSpecification_T, typename EvenPackInfo, typename OddPackInfo >
class CombinedInPlaceCpuPackInfo : public ::walberla::communication::UniformPackInfo
{
 public:
   template< typename... Args >
   CombinedInPlaceCpuPackInfo(std::shared_ptr< LatticeStorageSpecification_T >& storageSecification, Args&&... args)
      : storageSecification_(storageSecification), evenPackInfo_(std::forward< Args >(args)...), oddPackInfo_(std::forward< Args >(args)...)
   {}

   ~CombinedInPlaceCpuPackInfo() override = default;
   bool constantDataExchange() const override { return true; }
   bool threadsafeReceiving() const override { return true; }

   void unpackData(IBlock* receiver, stencil::Direction dir, mpi::RecvBuffer& buffer) override
   {
      if (storageSecification_->isEvenTimeStep())
      {
         return evenPackInfo_.unpackData(receiver, dir, buffer);
      }
      else
      {
         return oddPackInfo_.unpackData(receiver, dir, buffer);
      }
   }

   void communicateLocal(const IBlock* sender, IBlock* receiver, stencil::Direction dir) override
   {
      if (storageSecification_->isEvenTimeStep())
      {
         return evenPackInfo_.communicateLocal(sender, receiver, dir);
      }
      else
      {
         return oddPackInfo_.communicateLocal(sender, receiver, dir);
      }
   }

   void packDataImpl(const IBlock* sender, stencil::Direction dir, mpi::SendBuffer& outBuffer) const override
   {
      if (storageSecification_->isEvenTimeStep())
      {
         return evenPackInfo_.packDataImpl(sender, dir, outBuffer);
      }
      else
      {
         return oddPackInfo_.packDataImpl(sender, dir, outBuffer);
      }
   }

   void pack(stencil::Direction dir, unsigned char* buffer, IBlock* block) const
   {
      if (storageSecification_->isEvenTimeStep())
      {
         evenPackInfo_.pack(dir, buffer, block);
      }
      else
      {
         oddPackInfo_.pack(dir, buffer, block);
      }
   }

   void unpack(stencil::Direction dir, unsigned char* buffer, IBlock* block) const
   {
      if (storageSecification_->isEvenTimeStep())
      {
         evenPackInfo_.unpack(dir, buffer, block);
      }
      else
      {
         oddPackInfo_.unpack(dir, buffer, block);
      }
   }

   uint_t size(stencil::Direction dir, IBlock* block) const
   {
      if (storageSecification_->isEvenTimeStep())
      {
         return evenPackInfo_.size(dir, block);
      }
      else
      {
         return oddPackInfo_.size(dir, block);
      }
   }

 private:
   const std::shared_ptr< LatticeStorageSpecification_T >& storageSecification_;
   EvenPackInfo evenPackInfo_;
   OddPackInfo oddPackInfo_;
};

} // namespace walberla::lbm_generated
