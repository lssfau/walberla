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
//! \file PinnedMemoryBuffer.h
//! \ingroup cuda
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//! \brief Pinned Memory buffer for staging memory when using asynchronous CUDA memory copies.
//
//======================================================================================================================

#pragma once

#include "cuda/ErrorChecking.h"

#include <algorithm>
#include <cuda_runtime.h>


namespace walberla {
namespace cuda {
namespace communication {

template< typename T = unsigned char >
class GenericPinnedMemoryBuffer
{
public:
   typedef T ElementType;

   GenericPinnedMemoryBuffer();

   GenericPinnedMemoryBuffer( std::size_t initSize );

   GenericPinnedMemoryBuffer( const GenericPinnedMemoryBuffer & pb );

   ~GenericPinnedMemoryBuffer();

   inline T* ptr() const { return data_; }

   inline T* resize( std::size_t newSize );

   inline std::size_t size() const { return size_; }

   GenericPinnedMemoryBuffer & operator=( const GenericPinnedMemoryBuffer & pb ) = delete;
private:
   T * data_;
   std::size_t size_;
};

typedef GenericPinnedMemoryBuffer<> PinnedMemoryBuffer;

template< typename T >  // Element type
GenericPinnedMemoryBuffer<T>::GenericPinnedMemoryBuffer()
   : data_(nullptr), size_(0)
{
}


template< typename T >  // Element type
GenericPinnedMemoryBuffer<T>::GenericPinnedMemoryBuffer( std::size_t initSize )
   : data_(nullptr), size_(initSize)
{
   if (initSize > 0)
   {
      WALBERLA_CUDA_CHECK( cudaMallocHost( &data_, size_ * sizeof(T) ) );
   }
}

template< typename T >  // Element type
GenericPinnedMemoryBuffer<T>::GenericPinnedMemoryBuffer( const GenericPinnedMemoryBuffer & pb )
   : size_(pb.size_)
{
   if ( pb.size_ > 0 )
   {
      WALBERLA_CUDA_CHECK( cudaMallocHost( &data_, pb.size_ * sizeof(T) ) );

      std::copy( pb.data_, static_cast<T *>(pb.data_ + pb.size_), data_ );
   }
}

template< typename T >  // Element type
GenericPinnedMemoryBuffer<T>::~GenericPinnedMemoryBuffer()
{
   if ( data_ != nullptr )
   {
      WALBERLA_CUDA_CHECK( cudaFreeHost( data_ ) );
   }
}

template< typename T >  // Element type
T * GenericPinnedMemoryBuffer<T>::resize(std::size_t newSize)
{
   if ( newSize > size_ )
   {
      T * newBegin;
      WALBERLA_CUDA_CHECK( cudaMallocHost( &newBegin, newSize * sizeof(T) ) );

      std::swap( data_, newBegin );
      if ( newBegin != nullptr )
      {
         WALBERLA_CUDA_CHECK( cudaFreeHost( newBegin ) );
      }

      size_ = newSize;
   }

   return data_;
}

} // namespace communication
} // namespace cuda
} // namespace walberla
