
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
//! \file MPIDatatypes.impl.h
//! \ingroup field
//! \author Christian Godenschwager
//! \brief Implementation of functions to generate MPI data types for fields
//
//======================================================================================================================

#pragma once

namespace walberla {
namespace field {
namespace communication {


template<typename Field_T>
MPI_Datatype mpiDatatypeSlice( const Field_T & field,
                               const cell_idx_t xBeg, const cell_idx_t yBeg, const cell_idx_t zBeg, const cell_idx_t fBeg,
                               const cell_idx_t xEnd, const cell_idx_t yEnd, const cell_idx_t zEnd, const cell_idx_t fEnd )
{
   using T = typename Field_T::value_type;
   int sizes[4];
   int subsizes[4];
   int starts[4];

   if( field.layout() == field::fzyx )
   {
      sizes[0]    = int_c( field.fAllocSize() );
      sizes[1]    = int_c( field.zAllocSize() );
      sizes[2]    = int_c( field.yAllocSize() );
      sizes[3]    = int_c( field.xAllocSize() );

      subsizes[0] = int_c( fEnd - fBeg ) + 1;
      subsizes[1] = int_c( zEnd - zBeg ) + 1;
      subsizes[2] = int_c( yEnd - yBeg ) + 1;
      subsizes[3] = int_c( xEnd - xBeg ) + 1;

      starts[0]   = int_c( fBeg );
      starts[1]   = int_c( field.zOff() + zBeg );
      starts[2]   = int_c( field.yOff() + yBeg );
      starts[3]   = int_c( field.xOff() + xBeg );
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( field.layout(), field::zyxf );

      sizes[0]    = int_c( field.zAllocSize() );
      sizes[1]    = int_c( field.yAllocSize() );
      sizes[2]    = int_c( field.xAllocSize() );
      sizes[3]    = int_c( field.fAllocSize() );

      subsizes[0] = int_c( zEnd - zBeg ) + 1;
      subsizes[1] = int_c( yEnd - yBeg ) + 1;
      subsizes[2] = int_c( xEnd - xBeg ) + 1;
      subsizes[3] = int_c( fEnd - fBeg ) + 1;

      starts[0]   = int_c( field.zOff() + zBeg );
      starts[1]   = int_c( field.yOff() + yBeg );
      starts[2]   = int_c( field.xOff() + xBeg );
      starts[3]   = int_c( fBeg );
   }

   if( subsizes[0] < 1 || subsizes[1] < 1 || subsizes[2] < 1 || subsizes[3] < 1 )
   {
      // Create an empty dummy type for an empty interval
      MPI_Datatype emptyType;
      MPI_Type_contiguous( 0, MPITrait<T>::type(), &emptyType );
      return emptyType;
   }

   WALBERLA_DEBUG_SECTION()
   {
      for( int i = 0; i < 4; ++i )
      {
         WALBERLA_ASSERT_GREATER_EQUAL( subsizes[i], 1 );
         WALBERLA_ASSERT_LESS_EQUAL( subsizes[i], sizes[i] );
         WALBERLA_ASSERT_GREATER_EQUAL( starts[i], 0 );
         WALBERLA_ASSERT_LESS_EQUAL( starts[i], sizes[i] - subsizes[i] );
      }
   }

   MPI_Datatype newType;
   MPI_Type_create_subarray( 4, sizes, subsizes, starts, MPI_ORDER_C, MPITrait<T>::type(), &newType );

   return newType;
}



template<typename Field_T>
MPI_Datatype mpiDatatype( const Field_T & field )
{
   return mpiDatatypeSlice( field,
                            cell_idx_t( 0 ), cell_idx_t( 0 ), cell_idx_t( 0 ), cell_idx_t( 0 ),
                            cell_idx_t( field.xSize() ) - cell_idx_t( 1 ), cell_idx_t( field.ySize() ) - cell_idx_t( 1 ),
                            cell_idx_t( field.zSize() ) - cell_idx_t( 1 ), cell_idx_t( field.fSize() ) - cell_idx_t( 1 ) );
}


template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, cell_idx_t f /*= 0*/ )
{
   return mpiDatatypeSlice( field, 
                            interval.xMin(), interval.yMin(), interval.zMin(), f,
                            interval.xMax(), interval.yMax(), interval.zMax(), f );
}


template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, const cell_idx_t fBeg, const cell_idx_t fEnd )
{
   return mpiDatatypeSlice( field,
                            interval.xMin(), interval.yMin(), interval.zMin(), fBeg,
                            interval.xMax(), interval.yMax(), interval.zMax(), fEnd );
}


template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, const std::set<cell_idx_t> & fs )
{
   using T = typename Field_T::value_type;

   MPI_Datatype newType = MPI_DATATYPE_NULL;

   int sizes[3];
   int subsizes[3];
   int starts[3];

   sizes[0] = int_c( field.zAllocSize() );
   sizes[1] = int_c( field.yAllocSize() );
   sizes[2] = int_c( field.xAllocSize() );

   subsizes[0] = int_c( interval.zMax() - interval.zMin() ) + 1;
   subsizes[1] = int_c( interval.yMax() - interval.yMin() ) + 1;
   subsizes[2] = int_c( interval.xMax() - interval.xMin() ) + 1;

   starts[0] = int_c( field.zOff() + interval.zMin() );
   starts[1] = int_c( field.yOff() + interval.yMin() );
   starts[2] = int_c( field.xOff() + interval.xMin() );

   if( subsizes[0] < 1 || subsizes[1] < 1 || subsizes[2] < 1 )
   {
      // Create an empty dummy type for an empty interval
      MPI_Datatype emptyType;
      MPI_Type_contiguous( 0, MPITrait<T>::type(), &emptyType );
      return emptyType;
   }

   WALBERLA_DEBUG_SECTION()
   {
      for( int i = 0; i < 3; ++i )
      {
         WALBERLA_ASSERT_GREATER_EQUAL( starts[i], 0 );
         WALBERLA_ASSERT_LESS_EQUAL( starts[i], sizes[i] - subsizes[i] );
         WALBERLA_ASSERT_GREATER_EQUAL( subsizes[i], 1 );
         WALBERLA_ASSERT_LESS_EQUAL( subsizes[i], sizes[i] );
      }
   }

   if( field.layout() == field::fzyx )
   {
      MPI_Datatype tmpType = MPI_DATATYPE_NULL;
      MPI_Type_create_subarray( 3, sizes, subsizes, starts, MPI_ORDER_C, MPITrait<T>::type(), &tmpType );

      int count = int_c( fs.size() );
      std::vector<int> displacements( std::max( fs.size(), size_t( 1 ) ) ); // if "fs" is empty create a dummy vector from so that we can take an address to the first element
      std::transform( fs.begin(), fs.end(), displacements.begin(), int_c<cell_idx_t> );
      
      MPI_Type_create_indexed_block( count, 1, &( displacements.front() ), tmpType, &newType );
      
      MPI_Type_free( &tmpType );
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( field.layout(), field::zyxf );

      MPI_Datatype tmpType = MPI_DATATYPE_NULL;
      int count = int_c( fs.size() );
      std::vector<int> displacements( std::max( fs.size(), size_t(1) ) ); // if "fs" is empty create a dummy vector from so that we can take an address to the first element
      std::transform( fs.begin(), fs.end(), displacements.begin(), int_c<cell_idx_t> );

      MPI_Type_create_indexed_block( count, 1, &( displacements.front() ), MPITrait<T>::type(), &tmpType );

      MPI_Datatype resizedTmpType = MPI_DATATYPE_NULL;
      MPI_Type_create_resized( tmpType, 0, int_c( field.fAllocSize() * sizeof(T) ), &resizedTmpType );

      MPI_Type_create_subarray( 3, sizes, subsizes, starts, MPI_ORDER_C, resizedTmpType, &newType );

      MPI_Type_free( &tmpType );
      MPI_Type_free( &resizedTmpType );
   }

   return newType;
}


template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeWithGhostLayer( const GhostLayerField_T & field )
{
   return mpiDatatypeWithGhostLayer( field, field.nrOfGhostLayers() );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeWithGhostLayer( const GhostLayerField_T & field, const uint_t numGhostLayers )
{
   const cell_idx_t xBeg = - cell_idx_c( numGhostLayers );
   const cell_idx_t yBeg = - cell_idx_c( numGhostLayers );
   const cell_idx_t zBeg = - cell_idx_c( numGhostLayers );
   const cell_idx_t fBeg = cell_idx_t( 0 );

   const cell_idx_t xEnd = cell_idx_c( field.xSize() + numGhostLayers ) - cell_idx_t(1);
   const cell_idx_t yEnd = cell_idx_c( field.ySize() + numGhostLayers ) - cell_idx_t(1);
   const cell_idx_t zEnd = cell_idx_c( field.zSize() + numGhostLayers ) - cell_idx_t(1);
   const cell_idx_t fEnd = cell_idx_c( field.fSize() ) - cell_idx_t(1);

   return mpiDatatypeSlice( field,
                            xBeg, yBeg, zBeg, fBeg,
                            xEnd, yEnd, zEnd, fEnd );
}


template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnly( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice /*= false*/ )
{
   return mpiDatatypeGhostLayerOnly( field, field.nrOfGhostLayers(), dir, fullSlice );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnly( const GhostLayerField_T & field, const uint_t thickness, const stencil::Direction dir, const bool fullSlice /*= false*/ )
{
   CellInterval ci;
   field.getGhostRegion( dir, ci, cell_idx_c( thickness ), fullSlice );

   const cell_idx_t fBeg = cell_idx_t( 0 );
   const cell_idx_t fEnd = cell_idx_c( field.fSize() ) - cell_idx_t( 1 );

   return mpiDatatypeSliceXYZ( field, ci, fBeg, fEnd );
}


template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice /*= false*/, const cell_idx_t f /*= 0*/ )
{
   CellInterval ci;
   field.getGhostRegion( dir, ci, cell_idx_c( field.nrOfGhostLayers() ), fullSlice );

   return mpiDatatypeSliceXYZ( field, ci, f );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice, const cell_idx_t fBeg, const cell_idx_t fEnd )
{
   CellInterval ci;
   field.getGhostRegion( dir, ci, cell_idx_c( field.nrOfGhostLayers() ), fullSlice );

   return mpiDatatypeSliceXYZ( field, ci, fBeg, fEnd );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice, const std::set<cell_idx_t> & fs )
{
   CellInterval ci;
   field.getGhostRegion( dir, ci, cell_idx_c( field.nrOfGhostLayers() ), fullSlice );

   return mpiDatatypeSliceXYZ( field, ci, fs );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayer( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness /*= 1*/, const bool fullSlice /*= false*/ )
{
   CellInterval ci;
   field.getSliceBeforeGhostLayer( dir, ci, cell_idx_c( thickness ), fullSlice );

   const cell_idx_t fBeg = cell_idx_t( 0 );
   const cell_idx_t fEnd = cell_idx_c( field.fSize() ) - cell_idx_t( 1 );

   return mpiDatatypeSliceXYZ( field, ci, fBeg, fEnd );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness /*= 1*/, const cell_idx_t f /*= 0*/, const bool fullSlice /*= false*/ )
{
   CellInterval ci;
   field.getSliceBeforeGhostLayer( dir, ci, cell_idx_c( thickness ), fullSlice );

   return mpiDatatypeSliceXYZ( field, ci, f );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness, const cell_idx_t fBeg, const cell_idx_t fEnd, const bool fullSlice )
{
   CellInterval ci;
   field.getSliceBeforeGhostLayer( dir, ci, cell_idx_c( thickness ), fullSlice );

   return mpiDatatypeSliceXYZ( field, ci, fBeg, fEnd );
}

template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness, const std::set<cell_idx_t> & fs, const bool fullSlice )
{
   CellInterval ci;
   field.getSliceBeforeGhostLayer( dir, ci, cell_idx_c( thickness ), fullSlice );
   
   return mpiDatatypeSliceXYZ( field, ci, fs );
}


} // namespace communication
} // namespace field
} // namespace walberla
