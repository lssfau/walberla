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
//! \file BinaryRawFile.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"
#include "core/config/Config.h"

#include <fstream>
#include <string>
#include <vector>
#include <iterator>


namespace walberla {
namespace geometry {

class BinaryRawFile
{
public:
   BinaryRawFile( const std::string & filename, const Vector3< uint_t > & size, const std::string & datatype );

   template< typename T >
   static BinaryRawFile loadFile( const std::string & filename, const Vector3< uint_t > & size );

   BinaryRawFile( const Config::BlockHandle & configBlock );

   inline bool get( const Vector3< uint_t > & pos ) const;
   inline bool get( const uint_t x, const uint_t y, const uint_t z ) const;

   const Vector3< uint_t> & size() const { return size_; }

private:
   template< typename T >
   BinaryRawFile( const std::string & filename, const Vector3< uint_t > & size, const T );

   void init( const std::string & filename, const std::string & datatype );

   template< typename T >
   void init( const std::string & filename );

   Vector3< uint_t > size_;
   std::vector< bool > data_;
};

class BinaryRawFileInterpolator
{
public:
   enum Interpolator { NEAREST_NEIGHBOR };

   BinaryRawFileInterpolator( const AABB & aabb, const BinaryRawFile & binaryRawFile, Interpolator interpolator )
      : aabb_(aabb), interpolator_(interpolator), binaryRawFile_( binaryRawFile ) {}

   inline bool get( const Vector3< real_t > & pos ) const;
   inline bool get( const real_t x, const real_t y, const real_t z ) const;

   const AABB & aabb() const { return aabb_; }

   Interpolator interpolator() const { return interpolator_; }

private:
   inline bool getNearestNeighbor( const real_t x, const real_t y, const real_t z ) const;

   AABB aabb_;
   Interpolator interpolator_;
   const BinaryRawFile & binaryRawFile_;
};


template<typename T>
BinaryRawFile::BinaryRawFile( const std::string & filename, const Vector3<uint_t> & size, const T )
   : size_(size)
{
   init<T>( filename );
}


template<typename T>
void BinaryRawFile::init( const std::string & filename )
{
   const uint_t numElements = size_[0] * size_[1] * size_[2];
   data_.reserve( numElements );
   WALBERLA_ROOT_SECTION()
   {
      std::ifstream ifs( filename, std::ifstream::in | std::ifstream::binary );
      std::transform( std::istream_iterator<T>( ifs ), std::istream_iterator<T>(),
                      std::back_inserter( data_ ), []( const T v ) { return v > T( 0 ); } );
      WALBERLA_CHECK_EQUAL( numElements, data_.size(), "Error reading file \"" << filename << "\"!" );
   }
   mpi::broadcastObject( data_ );
}

template< typename T >
static BinaryRawFile loadFile( const std::string & filename, const Vector3< uint_t > & size )
{
   return BinaryRawFileReader( filename, size, T() );
}

bool BinaryRawFile::get( const Vector3< uint_t > & pos ) const
{
   return get( pos[0], pos[1], pos[2] );
}

inline bool BinaryRawFile::get( const uint_t x, const uint_t y, const uint_t z ) const
{
   WALBERLA_ASSERT_LESS( x, size_[0] );
   WALBERLA_ASSERT_LESS( y, size_[1] );
   WALBERLA_ASSERT_LESS( z, size_[2] );

   const uint_t i = z * size_[0] * size_[1] + y * size_[0] + x;
   
   WALBERLA_ASSERT_LESS( i, data_.size() );

   return data_[i];
}


bool BinaryRawFileInterpolator::get( const Vector3< real_t > & pos ) const
{
   return get( pos[0], pos[1], pos[2] );
}


bool BinaryRawFileInterpolator::get( const real_t x, const real_t y, const real_t z ) const
{
   WALBERLA_ASSERT( aabb_.contains( x, y, z ) );

   switch (interpolator_)
   {
   case NEAREST_NEIGHBOR:
      return getNearestNeighbor( x, y, z );
   default:
      WALBERLA_ABORT( "Unknown Interpolator!" );
   }
}

bool BinaryRawFileInterpolator::getNearestNeighbor( const real_t x, const real_t y, const real_t z ) const
{
   uint_t xInt = uint_c( (x - aabb_.xMin()) / aabb_.xSize() * real_t( binaryRawFile_.size()[0] ) + real_t( 0.5 ) );
   uint_t yInt = uint_c( (y - aabb_.yMin()) / aabb_.ySize() * real_t( binaryRawFile_.size()[1] ) + real_t( 0.5 ) );
   uint_t zInt = uint_c( (z - aabb_.zMin()) / aabb_.zSize() * real_t( binaryRawFile_.size()[2] ) + real_t( 0.5 ) );

   return binaryRawFile_.get( xInt, yInt, zInt );
}

} // namespace geometry
} // namespace walberla