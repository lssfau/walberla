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
//! \file SimpleCommunication.h
//! \ingroup blockforest
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/math/Vector3.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"

#include "lbm/communication/PdfFieldPackInfo.h"

namespace walberla
{
namespace blockforest
{
using communication::UniformBufferedScheme;

template< typename Stencil_T >
class SimpleCommunication : public communication::UniformBufferedScheme< Stencil_T >
{
   using RealScalarField_T      = GhostLayerField< real_t, 1 >;
   using VectorField_T          = GhostLayerField< Vector3< real_t >, 1 >;
   using VectorFieldFlattened_T = GhostLayerField< real_t, 3 >;
   using PdfField_T             = GhostLayerField< real_t, Stencil_T::Size >;
   using UintScalarField_T      = GhostLayerField< uint_t, 1 >;

   using FlagField16_T = FlagField< uint16_t >;
   using FlagField32_T = FlagField< uint32_t >;
   using FlagField64_T = FlagField< uint64_t >;

 public:
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3;
   }

   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3, BlockDataID f4)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3 << f4;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3, BlockDataID f4, BlockDataID f5)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3 << f4 << f5;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3, BlockDataID f4, BlockDataID f5, BlockDataID f6)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3 << f4 << f5 << f6;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3, BlockDataID f4, BlockDataID f5, BlockDataID f6, BlockDataID f7)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3 << f4 << f5 << f6 << f7;
   }
   SimpleCommunication(const std::weak_ptr< StructuredBlockForest >& blockForest, BlockDataID f1, BlockDataID f2,
                       BlockDataID f3, BlockDataID f4, BlockDataID f5, BlockDataID f6, BlockDataID f7, BlockDataID f8)
      : UniformBufferedScheme< Stencil_T >(blockForest), blockForest_(blockForest)
   {
      (*this) << f1 << f2 << f3 << f4 << f5 << f6 << f7 << f8;
   }

   SimpleCommunication& operator<<(BlockDataID fieldId)
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      if (blockForest->begin() == blockForest->end()) { return *this; }

      IBlock& firstBlock = *(blockForest->begin());

      using field::communication::PackInfo;

      if (firstBlock.isDataClassOrSubclassOf< PdfField_T >(fieldId))
      {
         this->addPackInfo(make_shared< PackInfo< PdfField_T > >(fieldId));
      }
      else
      {
         if (firstBlock.isDataClassOrSubclassOf< RealScalarField_T >(fieldId))
         {
            this->addPackInfo(make_shared< PackInfo< RealScalarField_T > >(fieldId));
         }
         else
         {
            if (firstBlock.isDataClassOrSubclassOf< VectorField_T >(fieldId))
            {
               this->addPackInfo(make_shared< PackInfo< VectorField_T > >(fieldId));
            }
            else
            {
               if (firstBlock.isDataClassOrSubclassOf< FlagField16_T >(fieldId))
               {
                  this->addPackInfo(make_shared< PackInfo< FlagField16_T > >(fieldId));
               }
               else
               {
                  if (firstBlock.isDataClassOrSubclassOf< FlagField32_T >(fieldId))
                  {
                     this->addPackInfo(make_shared< PackInfo< FlagField32_T > >(fieldId));
                  }
                  else
                  {
                     if (firstBlock.isDataClassOrSubclassOf< FlagField64_T >(fieldId))
                     {
                        this->addPackInfo(make_shared< PackInfo< FlagField64_T > >(fieldId));
                     }
                     else
                     {
                        if (firstBlock.isDataClassOrSubclassOf< UintScalarField_T >(fieldId))
                        {
                           this->addPackInfo(make_shared< PackInfo< UintScalarField_T > >(fieldId));
                        }
                        else
                        {
                           if (firstBlock.isDataClassOrSubclassOf< VectorFieldFlattened_T >(fieldId))
                           {
                              this->addPackInfo(make_shared< PackInfo< VectorFieldFlattened_T > >(fieldId));
                           }
                           else { WALBERLA_ABORT("Problem with UID"); }
                        }
                     }
                  }
               }
            }
         }
      }

      return *this;
   }

 protected:
   std::weak_ptr< StructuredBlockForest > blockForest_;
}; // class SimpleCommunication

} // namespace blockforest
} // namespace walberla
