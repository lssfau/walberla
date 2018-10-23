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
//! \file BasicBuffer.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


namespace walberla {
namespace cuda {
namespace communication {


   template<typename Allocator>
   CustomMemoryBuffer<Allocator>::CustomMemoryBuffer()
           : begin_( nullptr ), cur_( nullptr ), end_( nullptr )
   {
   }

   template<typename Allocator>
   CustomMemoryBuffer<Allocator>::CustomMemoryBuffer( std::size_t initSize )
           : begin_( nullptr ), cur_( nullptr ), end_( nullptr )
   {
      if( initSize > 0 )
      {
         begin_ = Allocator::allocate( initSize );
         end_ = begin_ + initSize;
         cur_ = begin_;
      }
   }

   template<typename Allocator>
   CustomMemoryBuffer<Allocator>::CustomMemoryBuffer( const CustomMemoryBuffer &pb )
   {
      if( pb.begin_ != nullptr )
      {
         begin_ = reinterpret_cast<ElementType *>(Allocator::allocate( pb.allocSize()));
         end_ = begin_ + pb.allocSize();
         Allocator::memcpy( begin_, pb.begin_, pb.allocSize());
         cur_ = begin_ + pb.size();
      }
   }

   template<typename Allocator>
   CustomMemoryBuffer<Allocator> &CustomMemoryBuffer<Allocator>::operator=( const CustomMemoryBuffer<Allocator> &pb )
   {
      auto copy( pb );
      std::swap( cur_, copy.cur_ );
      std::swap( begin_, copy.begin_ );
      std::swap( end_, copy.end_ );
      return *this;
   }

   template<typename Allocator>
   CustomMemoryBuffer<Allocator>::~CustomMemoryBuffer()
   {
      if( begin_ != nullptr )
         Allocator::deallocate( begin_ );
   }

   template<typename Allocator>
   void CustomMemoryBuffer<Allocator>::resize( std::size_t newSize )
   {
      if( newSize > allocSize())
      {
         auto offset = cur_ - begin_;

         ElementType *newBegin;

         newBegin = reinterpret_cast<ElementType *>(Allocator::allocate( newSize ));

         Allocator::memcpy( newBegin, begin_, size_t(end_ - begin_) );

         std::swap( begin_, newBegin );
         if( newBegin != nullptr )
            Allocator::deallocate( newBegin );

         end_ = begin_ + newSize;
         cur_ = begin_ + offset;
      }

   }

   template<typename Allocator>
   typename CustomMemoryBuffer<Allocator>::ElementType *CustomMemoryBuffer<Allocator>::advance( std::size_t bytes )
   {
      resize( size() + bytes );
      auto result = cur_;
      cur_ += bytes;
      WALBERLA_ASSERT_LESS_EQUAL( cur_, end_ );
      return result;
   }

   template<typename Allocator>
   typename CustomMemoryBuffer<Allocator>::ElementType *CustomMemoryBuffer<Allocator>::advanceNoResize( std::size_t bytes )
   {
      auto newSize = size() + bytes;
      if( newSize <= allocSize())
         return advance( bytes );
      else
         return nullptr;
   }

} // namespace communication
} // namespace cuda
} // namespace walberla
