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
//! \file DistributedSample.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <type_traits>
#include <iterator>
#include <vector>


namespace walberla {
namespace math {

/*******************************************************************************************************************//**
 * \brief   Class describing a distributed statistical sample
 **********************************************************************************************************************/
class DistributedSample
{
public:

   DistributedSample() :
      sum_( real_t(0) ), min_( real_t(0) ), max_( real_t(0) ), size_( uint_t(0) ), mean_( real_t(0) ), variance_( real_t(0) ) {}

   // insert

   void insert( const real_t val ) { data_.push_back( val ); }

   template <class InputIterator>
   void insert( InputIterator first, InputIterator last ) { for( auto it = first; it != last; ++it ) data_.push_back( *it ); }

   template< typename T >
   void castToRealAndInsert( const T & val );

   template <class InputIterator>
   void castToRealAndInsert( InputIterator first, InputIterator last );

   void clear() { data_.clear(); sum_ = real_t(0); min_ = real_t(0); max_ = real_t(0); size_ = uint_t(0); mean_ = real_t(0); variance_ = real_t(0); }

   // synchronization

   void mpiAllGather();
   void mpiGather( int rank );
   void mpiGatherRoot();

   // data retrieval

   real_t sum()   const { return sum_; }
   real_t min()   const { return min_; }
   real_t max()   const { return max_; }
   real_t range() const { return max_ - min_; }

   real_t mean() const { return mean_; }
   real_t avg() const { return mean(); }

   real_t variance() const { return variance_; }
   real_t stdDeviation() const { return std::sqrt( variance_ ); }
   real_t relativeStdDeviation() const { return std::sqrt( variance_ ) / mean_; }

   uint_t size() { return size_; }

   std::string format( const std::string & formatString = DEFAULT_FORMAT_STRING ) const;

private:
   
   std::vector< real_t > data_;

   real_t sum_;
   real_t min_;
   real_t max_;
   uint_t size_;
   real_t mean_;
   real_t variance_;

   static const std::string DEFAULT_FORMAT_STRING;
};


template< typename T >
void DistributedSample::castToRealAndInsert( const T & val )
{
   static_assert( std::is_arithmetic<T>::value, "You can only use DistributedSample::castToRealAndInsert with " \
                  "arithmetic types!" );

   insert( numeric_cast< real_t >( val ) );
}

template <class InputIterator>
void DistributedSample::castToRealAndInsert( InputIterator first, InputIterator last )
{
   static_assert( std::is_arithmetic< typename std::iterator_traits<InputIterator>::value_type >::value,
                  "You can only use DistributedSample::castToRealAndInsert with sequences of arithmetic types!" );

   while( first != last )
      castToRealAndInsert( *first++ );
}


} // namespace math
} // namespace walberla
