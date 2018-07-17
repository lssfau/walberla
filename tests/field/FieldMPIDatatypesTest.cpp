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
//! \file FieldMPIDatatypesTest.cpp
//! \ingroup field
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================


#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"

#include "field/Field.h"
#include "field/communication/MPIDatatypes.h"
#include "field/allocation/FieldAllocator.h"

#include "stencil/D3Q27.h"

#include <random>

namespace walberla {

using namespace field::communication;

class FieldRandomizer
{
public:

   template< typename T, uint_t fSize >
   void operator()( Field< T, fSize > & field )
   {
      std::uniform_real_distribution< T > distribution;

      for( auto it = field.begin(); it != field.end(); ++it )
      {
         *it = distribution( generator_ );
      }
   }

   template< typename T, uint_t fSize >
   void operator()( GhostLayerField< T, fSize > & field )
   {
      std::uniform_real_distribution< T > distribution;

      for( auto it = field.beginWithGhostLayer(); it != field.end(); ++it )
      {
         *it = distribution( generator_ );
      }
   }

private:
   std::mt19937 generator_;
};




template< typename SourceField, typename TargetField >
void testFullCopy( const SourceField & src, TargetField & dst )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   MPI_Datatype srcType = mpiDatatype( src );
   MPI_Datatype dstType = mpiDatatype( dst );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( src != dst );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( src == dst );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}

template< typename SourceField, typename TargetField >
void testOctantCopy( const SourceField & src, TargetField & dst )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   const CellInterval fullInterval = src.xyzSize();

   cell_idx_t xMid = fullInterval.xMin() + cell_idx_t( fullInterval.xSize() ) / cell_idx_t( 2 );
   cell_idx_t yMid = fullInterval.yMin() + cell_idx_t( fullInterval.ySize() ) / cell_idx_t( 2 );
   cell_idx_t zMid = fullInterval.zMin() + cell_idx_t( fullInterval.ySize() ) / cell_idx_t( 2 );

   CellInterval octants[] = {
      CellInterval( fullInterval.xMin(),  fullInterval.yMin(),  fullInterval.zMin(),  xMid,                 yMid,                 zMid                ),
      CellInterval( fullInterval.xMin(),  fullInterval.yMin(),  zMid + cell_idx_t(1), xMid,                 yMid,                 fullInterval.zMax() ),
      CellInterval( fullInterval.xMin(),  yMid + cell_idx_t(1), fullInterval.zMin(),  xMid,                 fullInterval.yMax(),  zMid                ),
      CellInterval( fullInterval.xMin(),  yMid + cell_idx_t(1), zMid + cell_idx_t(1), xMid,                 fullInterval.yMax(),  fullInterval.zMax() ),
      CellInterval( xMid + cell_idx_t(1), fullInterval.yMin(),  fullInterval.zMin(),  fullInterval.xMax(),  yMid,                 zMid                ),
      CellInterval( xMid + cell_idx_t(1), fullInterval.yMin(),  zMid + cell_idx_t(1), fullInterval.xMax(),  yMid,                 fullInterval.zMax() ),
      CellInterval( xMid + cell_idx_t(1), yMid + cell_idx_t(1), fullInterval.zMin(),  fullInterval.xMax(),  fullInterval.yMax(),  zMid                ),
      CellInterval( xMid + cell_idx_t(1), yMid + cell_idx_t(1), zMid + cell_idx_t(1), fullInterval.xMax(),  fullInterval.yMax(),  fullInterval.zMax() )
   };

   WALBERLA_CHECK( src != dst );

   for( auto it = octants; it != octants + 8; ++it )
   {
      MPI_Datatype srcType = mpiDatatypeSliceXYZ( src, *it, cell_idx_t( 0 ), cell_idx_c( src.fSize() ) - cell_idx_t( 1 ) );
      MPI_Datatype dstType = mpiDatatypeSliceXYZ( dst, *it, cell_idx_t( 0 ), cell_idx_c( dst.fSize() ) - cell_idx_t( 1 ) );

      MPI_Type_commit( &srcType );
      MPI_Type_commit( &dstType );
      
      MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

      MPI_Type_free( &srcType );
      MPI_Type_free( &dstType );
   }
   
   WALBERLA_CHECK( src == dst );
}


template< typename SourceField, typename TargetField >
void testSingleFCopy( const SourceField & src, TargetField & dst )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   WALBERLA_CHECK( src != dst );

   for( uint_t f = 0; f < src.fSize(); ++f )
   {
      MPI_Datatype srcType = mpiDatatypeSliceXYZ( src, src.xyzSize(), cell_idx_c( f ) );
      MPI_Datatype dstType = mpiDatatypeSliceXYZ( dst, dst.xyzSize(), cell_idx_c( f ) );

      MPI_Type_commit( &srcType );
      MPI_Type_commit( &dstType );

      MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

      MPI_Type_free( &srcType );
      MPI_Type_free( &dstType );
   }


   WALBERLA_CHECK( src == dst );
}


template< typename SourceField, typename TargetField >
void testSlicedFCopy( const SourceField & src, TargetField & dst )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   std::set<cell_idx_t> evenFs, oddFs;

   for( cell_idx_t f = 0; f < cell_idx_c( src.fSize() ); ++f )
   {
      if( f % cell_idx_t( 2 ) == 0 )
         evenFs.insert( f );
      else
         oddFs.insert( f );
   }

   WALBERLA_CHECK( src != dst );

   MPI_Datatype srcType = mpiDatatypeSliceXYZ( src, src.xyzSize(), evenFs );
   MPI_Datatype dstType = mpiDatatypeSliceXYZ( dst, dst.xyzSize(), evenFs );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );


   srcType = mpiDatatypeSliceXYZ( src, src.xyzSize(), oddFs );
   dstType = mpiDatatypeSliceXYZ( dst, dst.xyzSize(), oddFs );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );

   WALBERLA_CHECK( src == dst );
}


template< typename SourceField, typename TargetField >
void testIntervalCopy( const SourceField & src, TargetField & dst, const CellInterval & interval, const cell_idx_t fBeg, const cell_idx_t fEnd )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );
   
   dst.set( numeric_cast< typename SourceField::value_type >( 0 ) );

   WALBERLA_CHECK( src != dst );

   MPI_Datatype srcType = mpiDatatypeSliceXYZ( src, interval, fBeg, fEnd );
   MPI_Datatype dstType = mpiDatatypeSliceXYZ( dst, interval, fBeg, fEnd );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );

   CellInterval fieldInterval = src.xyzSize();
   for( auto it = fieldInterval.begin(); it != fieldInterval.end(); ++it )
   {
      if( interval.contains( *it ) )
      {
         for( cell_idx_t f = fBeg; f <= fEnd; ++f )
            WALBERLA_CHECK_IDENTICAL( src.get( *it, f ), dst.get( *it, f ) );
      }
      else
      {
         for( cell_idx_t f = fBeg; f <= fEnd; ++f )
            WALBERLA_CHECK_IDENTICAL( dst.get( *it, f ), numeric_cast< typename SourceField::value_type >( 0 ) );
      }
   }
}

template< typename SourceField, typename TargetField >
void testIntervalCopy( const SourceField & src, TargetField & dst )
{
   const CellInterval fieldInterval = src.xyzSize();

   for( cell_idx_t x = fieldInterval.xMin(); x <= fieldInterval.xMax(); ++x )
   {
      CellInterval testInterval = fieldInterval;
      testInterval.xMin() = testInterval.xMax() = x;

      testIntervalCopy( src, dst, testInterval, cell_idx_t( 0 ), cell_idx_c( src.fSize() ) - cell_idx_t( 1 ) );
   }

   for( cell_idx_t y = fieldInterval.yMin(); y <= fieldInterval.yMax(); ++y )
   {
      CellInterval testInterval = fieldInterval;
      testInterval.yMin() = testInterval.yMax() = y;

      testIntervalCopy( src, dst, testInterval, cell_idx_t( 0 ), cell_idx_c( src.fSize() ) - cell_idx_t( 1 ) );
   }

   for( cell_idx_t z = fieldInterval.zMin(); z <= fieldInterval.zMax(); ++z )
   {
      CellInterval testInterval = fieldInterval;
      testInterval.zMin() = testInterval.zMax() = z;

      testIntervalCopy( src, dst, testInterval, cell_idx_t( 0 ), cell_idx_c( src.fSize() ) - cell_idx_t( 1 ) );
   }

   for( cell_idx_t f = 0; f < cell_idx_c( src.fSize() ); ++f )
   {
      testIntervalCopy( src, dst, fieldInterval, f, f );
   }
}


template< typename SourceField, typename TargetField >
void runTests( SourceField & srcField, TargetField & dstField )
{
   dstField.set( 0.0 );
   testFullCopy( srcField, dstField );

   if( srcField.xSize() > uint_t( 1 ) && srcField.ySize() > uint_t( 1 ) && srcField.zSize() > uint_t( 1 ) )
   {
      dstField.set( 0.0 );
      testOctantCopy( srcField, dstField );
   }

   dstField.set( 0.0 );
   testSingleFCopy( srcField, dstField );

   dstField.set( 0.0 );
   testSlicedFCopy( srcField, dstField );

   dstField.set( 0.0 );
   testIntervalCopy( srcField, dstField );
   
}


template< typename SourceField, typename TargetField >
void testFullCopyGL( const SourceField & src, TargetField & dst )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   dst.setWithGhostLayer( 0.0 );

   MPI_Datatype srcType = mpiDatatypeWithGhostLayer( src );
   MPI_Datatype dstType = mpiDatatypeWithGhostLayer( dst );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( !std::equal( src.beginWithGhostLayer(), src.end(), dst.beginWithGhostLayer() ) );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( std::equal( src.beginWithGhostLayer(), src.end(), dst.beginWithGhostLayer() ) );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}

template< typename SourceField, typename TargetField >
void testPartialCopyGL( const SourceField & src, TargetField & dst, uint_t numberOfGhostLayers )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   dst.setWithGhostLayer( 0.0 );

   MPI_Datatype srcType = mpiDatatypeWithGhostLayer( src, numberOfGhostLayers );
   MPI_Datatype dstType = mpiDatatypeWithGhostLayer( dst, numberOfGhostLayers );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( !std::equal( src.beginWithGhostLayer( cell_idx_c( numberOfGhostLayers ) ), src.end(), dst.beginWithGhostLayer( cell_idx_c( numberOfGhostLayers ) ) ) );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( std::equal( src.beginWithGhostLayer( cell_idx_c( numberOfGhostLayers ) ), src.end(), dst.beginWithGhostLayer( cell_idx_c( numberOfGhostLayers ) ) ) );

   std::fill( dst.beginWithGhostLayer( cell_idx_c( numberOfGhostLayers ) ), dst.end(), typename TargetField::value_type( 0.0 ) );

   for( auto it = dst.beginWithGhostLayer(); it != dst.end(); ++it )
      WALBERLA_CHECK_IDENTICAL( *it, 0.0 );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}



template< typename SourceField, typename TargetField >
void testFullCopyGLOnly( const SourceField & src, TargetField & dst, const stencil::Direction dir, const bool fullSlice )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   dst.setWithGhostLayer( 0.0 );

   MPI_Datatype srcType = mpiDatatypeGhostLayerOnly( src, dir, fullSlice );
   MPI_Datatype dstType = mpiDatatypeGhostLayerOnly( dst, dir, fullSlice );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( !std::equal( src.beginGhostLayerOnly( dir, fullSlice ), src.end(), dst.beginGhostLayerOnly( dir, fullSlice ) ) );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( std::equal( src.beginGhostLayerOnly( dir, fullSlice ), src.end(), dst.beginGhostLayerOnly( dir, fullSlice ) ) );

   std::fill( dst.beginGhostLayerOnly( dir, fullSlice ), dst.end(), typename TargetField::value_type( 0.0 ) );

   for( auto it = dst.beginWithGhostLayer(); it != dst.end(); ++it )
      WALBERLA_CHECK_IDENTICAL( *it, 0.0 );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}

template< typename SourceField, typename TargetField >
void testPartialCopyGLOnly( const SourceField & src, TargetField & dst, const uint_t thickness, const stencil::Direction dir, const bool fullSlice )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   dst.setWithGhostLayer( 0.0 );

   MPI_Datatype srcType = mpiDatatypeGhostLayerOnly( src, thickness, dir, fullSlice );
   MPI_Datatype dstType = mpiDatatypeGhostLayerOnly( dst, thickness, dir, fullSlice );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( !std::equal( src.beginGhostLayerOnly( thickness, dir, fullSlice ), src.end(), dst.beginGhostLayerOnly( thickness, dir, fullSlice ) ) );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( std::equal( src.beginGhostLayerOnly( thickness, dir, fullSlice ), src.end(), dst.beginGhostLayerOnly( thickness, dir, fullSlice ) ) );

   std::fill( dst.beginGhostLayerOnly( thickness, dir, fullSlice ), dst.end(), typename TargetField::value_type( 0.0 ) );

   for( auto it = dst.beginWithGhostLayer(); it != dst.end(); ++it )
      WALBERLA_CHECK_IDENTICAL( *it, 0.0 );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}


template< typename SourceField, typename TargetField >
void testCopySliceBeforeGL( const SourceField & src, TargetField & dst, const uint_t thickness, const stencil::Direction dir, const bool fullSlice )
{
   WALBERLA_ASSERT( src.hasSameSize( dst ) );

   dst.setWithGhostLayer( 0.0 );

   MPI_Datatype srcType = mpiDatatypeSliceBeforeGhostlayer( src, dir, thickness, fullSlice );
   MPI_Datatype dstType = mpiDatatypeSliceBeforeGhostlayer( dst, dir, thickness, fullSlice );

   MPI_Type_commit( &srcType );
   MPI_Type_commit( &dstType );

   WALBERLA_CHECK( !std::equal( src.beginSliceBeforeGhostLayer( dir, cell_idx_c( thickness ), fullSlice ), src.end(), dst.beginSliceBeforeGhostLayer( dir, cell_idx_c( thickness ), fullSlice ) ) );

   MPI_Sendrecv( const_cast<typename SourceField::value_type *>( src.data() ), 1, srcType, 0, 0, dst.data(), 1, dstType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

   WALBERLA_CHECK( std::equal( src.beginSliceBeforeGhostLayer( dir, cell_idx_c( thickness ), fullSlice ), src.end(), dst.beginSliceBeforeGhostLayer( dir, cell_idx_c( thickness ), fullSlice ) ) );

   std::fill( dst.beginSliceBeforeGhostLayer( dir, cell_idx_c( thickness ), fullSlice ), dst.end(), typename TargetField::value_type( 0.0 ) );

   for( auto it = dst.beginWithGhostLayer(); it != dst.end(); ++it )
      WALBERLA_CHECK_IDENTICAL( *it, 0.0 );

   MPI_Type_free( &srcType );
   MPI_Type_free( &dstType );
}



template< typename SourceField, typename TargetField >
void runGhostLayerFieldTests( SourceField & srcField, TargetField & dstField )
{
   testFullCopyGL( srcField, dstField );  

   for( auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir )
   {
      testFullCopyGLOnly( srcField, dstField, *dir, true );
      testFullCopyGLOnly( srcField, dstField, *dir, false );
   }
   
   for( uint_t copyGl = 0; copyGl < srcField.nrOfGhostLayers(); ++copyGl )
   {
      testPartialCopyGL( srcField, dstField, copyGl );
   }

   for( uint_t copyGl = 1; copyGl < srcField.nrOfGhostLayers(); ++copyGl )
   {
      for( auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir )
      {
         testPartialCopyGLOnly( srcField, dstField, copyGl, *dir, true );
         testPartialCopyGLOnly( srcField, dstField, copyGl, *dir, false );
      }
   }

   for( auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir )
   {
      testCopySliceBeforeGL( srcField, dstField, 1, *dir, true );
      testCopySliceBeforeGL( srcField, dstField, 1, *dir, false );
   }
}

template< typename T, uint_t fSize >
void runTests( const Vector3<uint_t> & size, field::Layout layout, const shared_ptr< field::FieldAllocator<T> > & srcFieldAllocator, const shared_ptr< field::FieldAllocator<T> > & dstFieldAllocator )
{
   field::Field< T, fSize > srcField( size[0], size[1], size[2], 23.0, layout, srcFieldAllocator );
   field::Field< T, fSize > dstField( size[0], size[1], size[2], 42.0, layout, dstFieldAllocator );

   FieldRandomizer randomizer;
   randomizer( srcField );

   runTests( srcField, dstField );

   for( uint_t numGhostLayers = 1; numGhostLayers <= 3; ++numGhostLayers )
   {
      field::GhostLayerField< T, fSize > srcGhostLayerField( size[0], size[1], size[2], numGhostLayers, 23.0, layout, srcFieldAllocator );
      field::GhostLayerField< T, fSize > dstGhostLayerField( size[0], size[1], size[2], numGhostLayers, 42.0, layout, dstFieldAllocator );

      randomizer( srcGhostLayerField );

      runTests( srcGhostLayerField, dstGhostLayerField );
      runGhostLayerFieldTests( srcGhostLayerField, dstGhostLayerField );
   }
}

template< typename T, uint_t fSize >
void runTests( const Vector3<uint_t> & size )
{
   runTests<T, fSize>( size, field::fzyx, make_shared< field::StdFieldAlloc<T> >(),       make_shared< field::StdFieldAlloc<T>       >() );
   runTests<T, fSize>( size, field::zyxf, make_shared< field::StdFieldAlloc<T> >(),       make_shared< field::StdFieldAlloc<T>       >() );
   runTests<T, fSize>( size, field::fzyx, make_shared< field::AllocateAligned<T, 32> >(), make_shared< field::AllocateAligned<T, 32> >() );
   runTests<T, fSize>( size, field::zyxf, make_shared< field::AllocateAligned<T, 32> >(), make_shared< field::AllocateAligned<T, 32> >() );
   runTests<T, fSize>( size, field::fzyx, make_shared< field::StdFieldAlloc<T> >(),       make_shared< field::AllocateAligned<T, 32> >() );
   runTests<T, fSize>( size, field::zyxf, make_shared< field::StdFieldAlloc<T> >(),       make_shared< field::AllocateAligned<T, 32> >() );
   runTests<T, fSize>( size, field::fzyx, make_shared< field::AllocateAligned<T, 32> >(), make_shared< field::StdFieldAlloc<T>       >() );
   runTests<T, fSize>( size, field::zyxf, make_shared< field::AllocateAligned<T, 32> >(), make_shared< field::StdFieldAlloc<T>       >() );
}

template< typename T >
void runTests( const Vector3<uint_t> & size )
{
   runTests<T, 1 >( size );
   runTests<T, 2 >( size );
   runTests<T, 3 >( size );
}

void runTests( const Vector3<uint_t> & size )
{
   runTests< float  >( size );
   runTests< double >( size );
}

int main( int argc, char* argv[] )
{
   debug::enterTestMode();

   mpi::Environment mpiEnv( argc, argv );
   MPIManager::instance()->useWorldComm();
  
   runTests( Vector3<uint_t>( 1, 1, 1 ) );
   runTests( Vector3<uint_t>( 1, 1, 2 ) );
   runTests( Vector3<uint_t>( 1, 2, 1 ) );
   runTests( Vector3<uint_t>( 1, 2, 2 ) );
   runTests( Vector3<uint_t>( 2, 1, 1 ) );
   runTests( Vector3<uint_t>( 2, 1, 2 ) );
   runTests( Vector3<uint_t>( 2, 2, 1 ) );
   runTests( Vector3<uint_t>( 2, 2, 2 ) );

   runTests( Vector3<uint_t>( 1, 2, 3 ) );

   runTests( Vector3<uint_t>( 3, 3, 3 ) );
   runTests( Vector3<uint_t>( 4, 4, 4 ) );

   return EXIT_SUCCESS;
}
} // namespace walberla


int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
