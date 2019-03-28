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
//! \file Sample.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Declarations for class Sample
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/math/KahanSummation.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include <type_traits>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iterator>
#include <numeric>
#include <set>


namespace walberla {
namespace math {

/*******************************************************************************************************************//**
 * \brief   Class describing a statistical sample
 **********************************************************************************************************************/
class Sample : public std::multiset<real_t>
{
public:
   void merge( const Sample & other ) { insert( other.begin(), other.end() ); }

   template< typename T >
   iterator castToRealAndInsert(const T& val);

   template< typename T >
   iterator castToRealAndInsert(const_iterator position, const T& val);

   template <class InputIterator>
   void castToRealAndInsert(InputIterator first, InputIterator last);

   // collective MPI operations
   void mpiAllGather();
   void mpiGather(int rank);
   void mpiGatherRoot();

   // Data retrieval
   real_t sum()   const { WALBERLA_ASSERT(!empty()); return kahanSummation( begin(), end() ); }
   real_t min()   const { WALBERLA_ASSERT(!empty()); return *begin(); }
   real_t max()   const { WALBERLA_ASSERT(!empty()); return *--end(); }
   real_t range() const { return max() - min(); }

   real_t mean() const { WALBERLA_ASSERT(!empty()); return sum() / real_c(size()); }
   real_t median() const;

   real_t variance() const { return variance( mean() ); }
   real_t stdDeviation() const { return std::sqrt( variance() ); }
   real_t relativeStdDeviation() const;
   real_t mad() const;
   real_t giniCoefficient() const;

   real_t cummulativeDistributionFunction(const real_t x) const { WALBERLA_ASSERT(!empty()); return real_c(std::distance(begin(), upper_bound(x))) / real_c(size()); }
   real_t quantile(const real_t p) const;

   std::string format(const std::string & formatString = DEFAULT_FORMAT_STRING) const;

private:
   real_t variance( real_t mean ) const;
   
   static const std::string DEFAULT_FORMAT_STRING;
};

std::ostream & operator<<( std::ostream & os, const Sample & statReal );


template< typename T >
Sample::iterator Sample::castToRealAndInsert(const T& val)
{
   static_assert( std::is_arithmetic<T>::value, "You can only use Sample::castToRealAndInsert with " \
                  "arithmetic types!" );

   return insert( numeric_cast<value_type>( val ) );
}

template< typename T >
Sample::iterator Sample::castToRealAndInsert(const_iterator position, const T& val)
{
   static_assert( std::is_arithmetic<T>::value, "You can only use Sample::castToRealAndInsert with " \
                  "arithmetic types!" );

   return insert( position, numeric_cast<value_type>( val ) );
}

template <class InputIterator>
void Sample::castToRealAndInsert(InputIterator first, InputIterator last)
{
   static_assert( std::is_arithmetic< typename std::iterator_traits<InputIterator>::value_type >::value,
                  "You can only use Sample::castToRealAndInsert with sequences of arithmetic types!" );

   while( first != last )
      castToRealAndInsert( *first++ );
}


// Explicit call to operator<< due to problem with IBM compiler:
template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Sample & statReal )
{
   walberla::mpi::operator<<( buf, dynamic_cast< const std::multiset<real_t> & >( statReal ) );
   return buf;
}

template< typename T>    // Element type of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Sample & statReal )
{
   walberla::mpi::operator>>( buf, dynamic_cast< std::multiset<real_t> & >( statReal ) );
   return buf;
}

} // namespace math
} // namespace walberla
