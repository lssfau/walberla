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
//! \file GPUSweepBase.h
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================

#pragma once


#include "gpu/GPUField.h"

#include "core/debug/Debug.h"

#include "field/SwapableCompare.h"

#include <set>

namespace walberla {
namespace gpu
{


template < typename GPUField_T >
class GPUSweepBase
{
public:
   GPUSweepBase() = default;
   virtual ~GPUSweepBase()
   {
      for( auto field = dstFields_.begin(); field != dstFields_.end(); ++field )
      {
         delete *field;
      }
   }

   GPUField_T * getDstField( GPUField_T * const src )
   {
      auto it = dstFields_.find( src );
      if( it != dstFields_.end() )
      {
         return *it;
      }

      GPUField_T * dst = src->cloneUninitialized();
      WALBERLA_ASSERT_NOT_NULLPTR( dst )

      dstFields_.insert( dst );

      return dst;
   }

protected:

   std::set< GPUField_T *, field::SwapableCompare< GPUField_T * > > dstFields_;
};


} // namespace gpu
} // namespace walberla

