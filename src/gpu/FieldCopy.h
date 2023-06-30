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
//! \file FieldCopy.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "ErrorChecking.h"
#include "GPUField.h"

#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/Field.h"
#include "field/GhostLayerField.h"

#include "core/Abort.h"
#include "core/logging/Logging.h"


namespace walberla {
namespace gpu
{


   template<typename DstType, typename SrcType>
   void fieldCpy( const shared_ptr< StructuredBlockStorage > & blocks,  BlockDataID dstID, ConstBlockDataID srcID )
   {
      for ( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
               DstType * dst = blockIt->getData<DstType>( dstID );
         const SrcType * src = blockIt->getData<SrcType>( srcID );
         fieldCpy( *dst, *src );
      }
   }

   template<typename DstType, typename SrcType>
   std::function<void()> fieldCpyFunctor( const shared_ptr< StructuredBlockStorage > & blocks,
                                            BlockDataID dstID, ConstBlockDataID srcID )
   {
      return std::bind( fieldCpy<DstType,SrcType>, blocks, dstID, srcID );
   }



   template<typename DstType, typename SrcType>
   void fieldCpySweepFunction( BlockDataID dstID, ConstBlockDataID srcID, IBlock * block )
   {
            DstType * dst = block->getData<DstType>( dstID );
      const SrcType * src = block->getData<SrcType>( srcID );
      fieldCpy( *dst, *src );
   }

   template<typename DstType, typename SrcType>
   std::function<void(IBlock*)> fieldCpyFunctor( BlockDataID dstID, ConstBlockDataID srcID )
   {
      return std::bind( fieldCpySweepFunction<DstType,SrcType>, dstID, srcID, std::placeholders::_1 );
   }





   template<typename T, uint_t fs>
   void fieldCpy(gpu::GPUField<T> & dst, const field::Field<T,fs> & src );



   template<typename T, uint_t fs>
   void fieldCpy( field::Field<T,fs> & dst, const gpu::GPUField<T> & src );




   //===================================================================================================================
   //
   //  Implementation
   //
   //===================================================================================================================




   template<typename T, uint_t fs>
   void fieldCpy(gpu::GPUField<T> & dst, const field::Field<T,fs> & src )
   {
      WALBERLA_DEVICE_SECTION()
      {
         gpuMemcpy3DParms p;
         memset(&p, 0, sizeof(p));

         if (dst.layout() != src.layout()) { WALBERLA_ABORT("Cannot copy fields with different layout") }

         bool canCopy =
            (src.layout() == fzyx && dst.fAllocSize() == src.fAllocSize() && dst.zAllocSize() == src.zAllocSize() &&
             dst.yAllocSize() == src.yAllocSize() && dst.xSize() == src.xSize()) ||
            (src.layout() == zyxf && dst.zAllocSize() == src.zAllocSize() && dst.yAllocSize() == src.yAllocSize() &&
             dst.xAllocSize() == src.xAllocSize() && dst.fSize() == src.fSize());

         if (!canCopy) { WALBERLA_ABORT("Field have to have the same size ") }

         if (dst.layout() == fzyx)
         {
            p.srcPtr = make_gpuPitchedPtr((void*) (src.data()),         // pointer
                                          sizeof(T) * src.xAllocSize(), // pitch
                                          src.xAllocSize(),             // inner dimension size
                                          src.yAllocSize());            // next outer dimension size

            p.extent.width  = std::min(dst.xAllocSize(), src.xAllocSize()) * sizeof(T);
            p.extent.height = dst.yAllocSize();
            p.extent.depth  = dst.zAllocSize() * dst.fAllocSize();
         }
         else
         {
            p.srcPtr = make_gpuPitchedPtr((void*) (src.data()),         // pointer
                                          sizeof(T) * src.fAllocSize(), // pitch
                                          src.fAllocSize(),             // inner dimension size
                                          src.xAllocSize());            // next outer dimension size

            p.extent.width  = std::min(dst.fAllocSize(), src.fAllocSize()) * sizeof(T);
            p.extent.height = dst.xAllocSize();
            p.extent.depth  = dst.yAllocSize() * dst.zAllocSize();
         }

         p.dstPtr = dst.pitchedPtr();
         p.kind   = gpuMemcpyHostToDevice;
         WALBERLA_GPU_CHECK(gpuMemcpy3D(&p))
      }
   }



   template<typename T, uint_t fs>
   void fieldCpy( field::Field<T,fs> & dst, const gpu::GPUField<T> & src )
   {
      WALBERLA_DEVICE_SECTION()
      {
         gpuMemcpy3DParms p;
         memset(&p, 0, sizeof(p));

         if (dst.layout() != src.layout()) { WALBERLA_ABORT("Cannot copy fields with different layout") }

         bool canCopy =
            (src.layout() == fzyx && dst.fAllocSize() == src.fAllocSize() && dst.zAllocSize() == src.zAllocSize() &&
             dst.yAllocSize() == src.yAllocSize() && dst.xSize() == src.xSize()) ||
            (src.layout() == zyxf && dst.zAllocSize() == src.zAllocSize() && dst.yAllocSize() == src.yAllocSize() &&
             dst.xAllocSize() == src.xAllocSize() && dst.fSize() == src.fSize());

         if (!canCopy) { WALBERLA_ABORT("Field have to have the same size ") }

         if (dst.layout() == fzyx)
         {
            p.dstPtr = make_gpuPitchedPtr((void*) (dst.data()),         // pointer
                                          sizeof(T) * dst.xAllocSize(), // pitch
                                          dst.xAllocSize(),             // inner dimension size
                                          dst.yAllocSize());            // next outer dimension size

            p.extent.width  = std::min(dst.xAllocSize(), src.xAllocSize()) * sizeof(T);
            p.extent.height = dst.yAllocSize();
            p.extent.depth  = dst.zAllocSize() * dst.fAllocSize();
         }
         else
         {
            p.dstPtr = make_gpuPitchedPtr((void*) (dst.data()),         // pointer
                                          sizeof(T) * dst.fAllocSize(), // pitch
                                          dst.fAllocSize(),             // inner dimension size
                                          dst.xAllocSize());            // next outer dimension size

            p.extent.width  = std::min(dst.fAllocSize(), src.fAllocSize()) * sizeof(T);
            p.extent.height = dst.xAllocSize();
            p.extent.depth  = dst.yAllocSize() * dst.zAllocSize();
         }

         p.srcPtr = src.pitchedPtr();
         p.kind   = gpuMemcpyDeviceToHost;
         WALBERLA_GPU_CHECK(gpuMemcpy3D(&p))
      }
   }

} // namespace gpu
} // namespace walberla


