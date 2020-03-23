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
//! \file BlockDataHandling.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockDataHandling.h"
#include "blockforest/StructuredBlockForest.h"
#include "core/debug/CheckFunctions.h"
#include "core/math/Vector2.h"
#include "core/math/Vector3.h"
#include "field/FlagField.h"

#include <type_traits>


namespace walberla {
namespace field {



// still virtual, one must implement protected member functions 'allocate' and 'reallocate'
template< typename Field_T, bool Pseudo2D = false >
class BlockDataHandling : public blockforest::BlockDataHandling< Field_T >
{
public:

   typedef typename Field_T::value_type Value_T;
   typedef std::function< void ( Field_T * field, IBlock * const block ) > InitializationFunction_T;

   virtual ~BlockDataHandling() {}

   void addInitializationFunction( const InitializationFunction_T & initFunction ) { initFunction_ = initFunction; }

   Field_T * initialize( IBlock * const block )
   {
      Field_T * field = allocate( block );
      
      if( initFunction_ )
         initFunction_( field, block );

      return field;
   }

   inline void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer );

   void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child );
   void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer );

   Field_T * deserialize( IBlock * const block ) { return reallocate( block ); }

   Field_T * deserializeCoarseToFine( Block * const block ) { return reallocate( block ); }
   Field_T * deserializeFineToCoarse( Block * const block ) { return reallocate( block ); }   
   
   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer );

   void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer );
   void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child );

protected:

   /// must be thread-safe !
   virtual Field_T *   allocate( IBlock * const block ) = 0; // used in 'initialize'
   /// must be thread-safe !
   virtual Field_T * reallocate( IBlock * const block ) = 0; // used in all deserialize member functions

   template< typename T > struct Merge
   { static T result( const T & value ) { return Pseudo2D ? static_cast<T>( value / numeric_cast<T>(4) ) : static_cast<T>( value / numeric_cast<T>(8) ); } };

   template< typename T > struct Merge< Vector2<T> >
   { static Vector2<T> result( const Vector2<T> & value ) { return Pseudo2D ? (value / numeric_cast<T>(4)) : (value / numeric_cast<T>(8)); } };

   template< typename T > struct Merge< Vector3<T> >
   { static Vector3<T> result( const Vector3<T> & value ) { return Pseudo2D ? (value / numeric_cast<T>(4)) : (value / numeric_cast<T>(8)); } };

   void sizeCheck( const uint_t xSize, const uint_t ySize, const uint_t zSize )
   {
      WALBERLA_CHECK( (xSize & uint_t(1)) == uint_t(0), "The x-size of your field must be divisible by 2." );
      WALBERLA_CHECK( (ySize & uint_t(1)) == uint_t(0), "The y-size of your field must be divisible by 2." );
      if( Pseudo2D )
      { WALBERLA_CHECK( zSize == uint_t(1), "The z-size of your field must be equal to 1 (pseudo 2D mode)." ); }
      else
      { WALBERLA_CHECK( (zSize & uint_t(1)) == uint_t(0), "The z-size of your field must be divisible by 2." ); }
   }
   
   InitializationFunction_T initFunction_;

}; // class BlockDataHandling



template< typename Field_T, bool Pseudo2D >
inline void BlockDataHandling< Field_T, Pseudo2D >::serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
{
   Field_T * field = block->template getData< Field_T >(id);
   WALBERLA_ASSERT_NOT_NULLPTR( field );

#ifndef NDEBUG
   buffer << field->xSize() << field->ySize() << field->zSize() << field->fSize();
#endif

   for( auto it = field->begin(); it != field->end(); ++it )
      buffer << *it;
}



template< typename Field_T, bool Pseudo2D >
void BlockDataHandling< Field_T, Pseudo2D >::serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child )
{
   Field_T * field = block->template getData< Field_T >(id);
   WALBERLA_ASSERT_NOT_NULLPTR( field );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   const uint_t fSize = field->fSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   buffer << child << ( xSize / uint_t(2) ) << ( ySize / uint_t(2) ) << ( Pseudo2D ? zSize : ( zSize / uint_t(2) ) ) << fSize;
#endif

   const cell_idx_t zBegin = Pseudo2D ? cell_idx_t(0) : ( (child & uint_t(4)) ? ( cell_idx_c( zSize ) / cell_idx_t(2) ) : cell_idx_t(0) );
   const cell_idx_t zEnd = Pseudo2D ? cell_idx_t(1) : ( (child & uint_t(4)) ? cell_idx_c( zSize ) : ( cell_idx_c( zSize ) / cell_idx_t(2) ) );
   for( cell_idx_t z = zBegin; z < zEnd; ++z )
   {
      const cell_idx_t yEnd = (child & uint_t(2)) ? cell_idx_c( ySize ) : ( cell_idx_c( ySize ) / cell_idx_t(2) );
      for( cell_idx_t y = (child & uint_t(2)) ? ( cell_idx_c( ySize ) / cell_idx_t(2) ) : cell_idx_t(0); y < yEnd; ++y )
      {
         const cell_idx_t xEnd = (child & uint_t(1)) ? cell_idx_c( xSize ) : ( cell_idx_c( xSize ) / cell_idx_t(2) );
         for( cell_idx_t x = (child & uint_t(1)) ? ( cell_idx_c( xSize ) / cell_idx_t(2) ) : cell_idx_t(0); x < xEnd; ++x )
         {
            for( uint_t f = uint_t(0); f < fSize; ++f )
               buffer << field->get(x,y,z,f);
         }
      }
   }
}



template< typename Field_T, bool Pseudo2D >
void BlockDataHandling< Field_T, Pseudo2D >::serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
{
   Field_T * field = block->template getData< Field_T >(id);
   WALBERLA_ASSERT_NOT_NULLPTR( field );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   const uint_t fSize = field->fSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   buffer << block->getId().getBranchId() << ( xSize / uint_t(2) ) << ( ySize / uint_t(2) ) << ( Pseudo2D ? zSize : ( zSize / uint_t(2) ) ) << fSize;
#endif

   for( cell_idx_t z = cell_idx_t(0); z < cell_idx_c( zSize ); z += cell_idx_t(2) ) {
      for( cell_idx_t y = cell_idx_t(0); y < cell_idx_c( ySize ); y += cell_idx_t(2) ) {
         for( cell_idx_t x = cell_idx_t(0); x < cell_idx_c( xSize ); x += cell_idx_t(2) ) {
            for( uint_t f = uint_t(0); f < fSize; ++f )
            {
               Value_T result =                                  field->get( x,                 y,                 z,                 f );
                       result = static_cast< Value_T >( result + field->get( x + cell_idx_t(1), y                , z                , f ) );
                       result = static_cast< Value_T >( result + field->get( x                , y + cell_idx_t(1), z                , f ) );
                       result = static_cast< Value_T >( result + field->get( x + cell_idx_t(1), y + cell_idx_t(1), z                , f ) );
               if( ! Pseudo2D )
               {
                       result = static_cast< Value_T >( result + field->get( x                , y                , z + cell_idx_t(1), f ) );
                       result = static_cast< Value_T >( result + field->get( x + cell_idx_t(1), y                , z + cell_idx_t(1), f ) );
                       result = static_cast< Value_T >( result + field->get( x                , y + cell_idx_t(1), z + cell_idx_t(1), f ) );
                       result = static_cast< Value_T >( result + field->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1), f ) );
               }

               buffer << Merge< Value_T >::result( result );
            }
         }
      }
   }
}



template< typename Field_T, bool Pseudo2D >
inline void BlockDataHandling< Field_T, Pseudo2D >::deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   Field_T * field = block->template getData< Field_T >( id );

#ifndef NDEBUG
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   uint_t fSender( uint_t(0) );
   buffer >> xSender >> ySender >> zSender >> fSender;
   WALBERLA_ASSERT_EQUAL( xSender, field->xSize() );
   WALBERLA_ASSERT_EQUAL( ySender, field->ySize() );
   WALBERLA_ASSERT_EQUAL( zSender, field->zSize() );
   WALBERLA_ASSERT_EQUAL( fSender, field->fSize() );
#endif

   for( auto it = field->begin(); it != field->end(); ++it )
      buffer >> *it;
}



template< typename Field_T, bool Pseudo2D >
void BlockDataHandling< Field_T, Pseudo2D >::deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   Field_T * field = block->template getData< Field_T >( id );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   const uint_t fSize = field->fSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   uint_t branchId( uint_t(0) );
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   uint_t fSender( uint_t(0) );
   buffer >> branchId >> xSender >> ySender >> zSender >> fSender;
   WALBERLA_ASSERT_EQUAL( branchId, block->getId().getBranchId() );
   WALBERLA_ASSERT_EQUAL( xSender, xSize / uint_t(2) );
   WALBERLA_ASSERT_EQUAL( ySender, ySize / uint_t(2) );
   if( Pseudo2D )
   { WALBERLA_ASSERT_EQUAL( zSender, zSize ); }
   else
   { WALBERLA_ASSERT_EQUAL( zSender, zSize / uint_t(2) ); }
   WALBERLA_ASSERT_EQUAL( fSender, fSize );
#endif

   for( cell_idx_t z = cell_idx_t(0); z < cell_idx_c( zSize ); z += cell_idx_t(2) ) {
      for( cell_idx_t y = cell_idx_t(0); y < cell_idx_c( ySize ); y += cell_idx_t(2) ) {
         for( cell_idx_t x = cell_idx_t(0); x < cell_idx_c( xSize ); x += cell_idx_t(2) ) {
            for( uint_t f = uint_t(0); f < fSize; ++f )
            {
               Value_T value;
               buffer >> value;

               field->get( x,                 y,                 z,                 f ) = value;
               field->get( x + cell_idx_t(1), y                , z                , f ) = value;
               field->get( x                , y + cell_idx_t(1), z                , f ) = value;
               field->get( x + cell_idx_t(1), y + cell_idx_t(1), z                , f ) = value;
               if( ! Pseudo2D )
               {
                  field->get( x                , y                , z + cell_idx_t(1), f ) = value;
                  field->get( x + cell_idx_t(1), y                , z + cell_idx_t(1), f ) = value;
                  field->get( x                , y + cell_idx_t(1), z + cell_idx_t(1), f ) = value;
                  field->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1), f ) = value;
               }
            }
         }
      }
   }
}



template< typename Field_T, bool Pseudo2D >
void BlockDataHandling< Field_T, Pseudo2D >::deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child )
{
   Field_T * field = block->template getData< Field_T >( id );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   const uint_t fSize = field->fSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   uint_t branchId( uint_t(0) );
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   uint_t fSender( uint_t(0) );
   buffer >> branchId >> xSender >> ySender >> zSender >> fSender;
   WALBERLA_ASSERT_EQUAL( branchId, child );
   WALBERLA_ASSERT_EQUAL( xSender, xSize / uint_t(2) );
   WALBERLA_ASSERT_EQUAL( ySender, ySize / uint_t(2) );
   if( Pseudo2D )
   { WALBERLA_ASSERT_EQUAL( zSender, zSize ); }
   else
   { WALBERLA_ASSERT_EQUAL( zSender, zSize / uint_t(2) ); }
   WALBERLA_ASSERT_EQUAL( fSender, fSize );
#endif

   const cell_idx_t zBegin = Pseudo2D ? cell_idx_t(0) : ( (child & uint_t(4)) ? ( cell_idx_c( zSize ) / cell_idx_t(2) ) : cell_idx_t(0) );
   const cell_idx_t zEnd = Pseudo2D ? cell_idx_t(1) : ( (child & uint_t(4)) ? cell_idx_c( zSize ) : ( cell_idx_c( zSize ) / cell_idx_t(2) ) );
   for( cell_idx_t z = zBegin; z < zEnd; ++z )
   {
      const cell_idx_t yEnd = (child & uint_t(2)) ? cell_idx_c( ySize ) : ( cell_idx_c( ySize ) / cell_idx_t(2) );
      for( cell_idx_t y = (child & uint_t(2)) ? ( cell_idx_c( ySize ) / cell_idx_t(2) ) : cell_idx_t(0); y < yEnd; ++y )
      {
         const cell_idx_t xEnd = (child & uint_t(1)) ? cell_idx_c( xSize ) : ( cell_idx_c( xSize ) / cell_idx_t(2) );
         for( cell_idx_t x = (child & uint_t(1)) ? ( cell_idx_c( xSize ) / cell_idx_t(2) ) : cell_idx_t(0); x < xEnd; ++x )
         {
            for( uint_t f = uint_t(0); f < fSize; ++f )
               buffer >> field->get(x,y,z,f);
         }
      }
   }
}






// allocation helper functions used in class 'DefaultBlockDataHandling' (see below)
namespace internal
{

template< typename GhostLayerField_T >
inline GhostLayerField_T * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl,
                                     const typename GhostLayerField_T::value_type & v, Layout l )
{
   return new GhostLayerField_T(x,y,z,gl,v,l);
}
template<>
inline FlagField<uint8_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, const uint8_t &, Layout )
{
   return new FlagField<uint8_t>(x,y,z,gl);
}
template<>
inline FlagField<uint16_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, const uint16_t &, Layout )
{
   return new FlagField<uint16_t>(x,y,z,gl);
}
template<>
inline FlagField<uint32_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, const uint32_t &, Layout )
{
   return new FlagField<uint32_t>(x,y,z,gl);
}
template<>
inline FlagField<uint64_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, const uint64_t &, Layout )
{
   return new FlagField<uint64_t>(x,y,z,gl);
}

template< typename GhostLayerField_T >
inline GhostLayerField_T * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, Layout l )
{
   return new GhostLayerField_T(x,y,z,gl,l);
}
template<>
inline FlagField<uint8_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, Layout )
{
   return new FlagField<uint8_t>(x,y,z,gl);
}
template<>
inline FlagField<uint16_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, Layout )
{
   return new FlagField<uint16_t>(x,y,z,gl);
}
template<>
inline FlagField<uint32_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, Layout )
{
   return new FlagField<uint32_t>(x,y,z,gl);
}
template<>
inline FlagField<uint64_t> * allocate( const uint_t x, const uint_t y, const uint_t z, const uint_t gl, Layout )
{
   return new FlagField<uint64_t>(x,y,z,gl);
}

inline Vector3< uint_t > defaultSize( const shared_ptr< StructuredBlockStorage > & blocks, IBlock * const block )
{
   return Vector3<uint_t>( blocks->getNumberOfXCells( *block ), blocks->getNumberOfYCells( *block ), blocks->getNumberOfZCells( *block ) );
}

} // namespace internal



template< typename GhostLayerField_T >
class DefaultBlockDataHandling : public BlockDataHandling< GhostLayerField_T >
{
public:

   typedef typename GhostLayerField_T::value_type Value_T;

   DefaultBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks,
                             const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( uint_t(1) ), initValue_(), layout_( zyxf ), calculateSize_( calculateSize )
   {}

   DefaultBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks, const uint_t nrOfGhostLayers,
                             const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( nrOfGhostLayers ), initValue_(), layout_( zyxf ), calculateSize_( calculateSize )
   {}

   DefaultBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks, const uint_t nrOfGhostLayers,
                             const Value_T & initValue, const Layout layout = zyxf,
                             const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( nrOfGhostLayers ), initValue_( initValue ), layout_( layout ), calculateSize_( calculateSize )
   {
      static_assert( !std::is_same< GhostLayerField_T, FlagField< Value_T > >::value,
                     "When using class FlagField, only constructors without the explicit specification of an initial value and the field layout are available!" );
   }

protected:

   GhostLayerField_T * allocate( IBlock * const block )
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'DefaultBlockDataHandling' for a block storage object that doesn't exist anymore" );
      const Vector3< uint_t > size = calculateSize_( blocks, block );
      return internal::allocate< GhostLayerField_T >( size[0], size[1], size[2],
                                                      nrOfGhostLayers_, initValue_, layout_ );
   }

   GhostLayerField_T * reallocate( IBlock * const block )
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'DefaultBlockDataHandling' for a block storage object that doesn't exist anymore" );
      const Vector3< uint_t > size = calculateSize_( blocks, block );
      return internal::allocate< GhostLayerField_T >( size[0], size[1], size[2],
                                                      nrOfGhostLayers_, layout_ );
   }

private:

   weak_ptr< StructuredBlockStorage > blocks_;

   uint_t  nrOfGhostLayers_;
   Value_T initValue_;
   Layout  layout_;
   const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize_;

}; // class DefaultBlockDataHandling






template< typename GhostLayerField_T >
class AlwaysInitializeBlockDataHandling : public blockforest::AlwaysInitializeBlockDataHandling< GhostLayerField_T >
{
public:

   typedef typename GhostLayerField_T::value_type Value_T;
   typedef std::function< void ( GhostLayerField_T * field, IBlock * const block ) > InitializationFunction_T;

   AlwaysInitializeBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks,
                                      const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( uint_t(1) ), initValue_(), layout_( zyxf ), calculateSize_( calculateSize )
   {}

   AlwaysInitializeBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks, const uint_t nrOfGhostLayers,
                                      const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( nrOfGhostLayers ), initValue_(), layout_( zyxf ), calculateSize_( calculateSize )
   {}

   AlwaysInitializeBlockDataHandling( const weak_ptr< StructuredBlockStorage > & blocks, const uint_t nrOfGhostLayers,
                                      const Value_T & initValue, const Layout layout,
                                      const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = internal::defaultSize ) :
      blocks_( blocks ), nrOfGhostLayers_( nrOfGhostLayers ), initValue_( initValue ), layout_( layout ), calculateSize_( calculateSize )
   {
      static_assert( ! std::is_same< GhostLayerField_T, FlagField< Value_T > >::value,
                     "When using class FlagField, only constructors without the explicit specification of an initial value and the field layout are available!" );
   }

   void addInitializationFunction( const InitializationFunction_T & initFunction ) { initFunction_ = initFunction; }

   GhostLayerField_T * initialize( IBlock * const block )
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'AlwaysInitializeBlockDataHandling' for a block storage object that doesn't exist anymore" );
      Vector3<uint_t> size = calculateSize_( blocks, block );
      GhostLayerField_T * field = internal::allocate< GhostLayerField_T >( size[0], size[1], size[2],
                                                                           nrOfGhostLayers_, initValue_, layout_ );
      if( initFunction_ )
         initFunction_( field, block );

      return field;
   }

private:

   weak_ptr< StructuredBlockStorage > blocks_;

   uint_t  nrOfGhostLayers_;
   Value_T initValue_;
   Layout  layout_;
   const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize_;

   InitializationFunction_T initFunction_;

}; // class AlwaysInitializeBlockDataHandling






template< typename Field_T >
class CloneBlockDataHandling : public blockforest::AlwaysInitializeBlockDataHandling< Field_T >
{
public:

   CloneBlockDataHandling( const ConstBlockDataID & fieldToClone ) :
      fieldToClone_( fieldToClone )
   {}

   Field_T * initialize( IBlock * const block )
   {
      const Field_T * toClone = block->template getData< Field_T >( fieldToClone_ );
      return toClone->clone();
   }

private:

   ConstBlockDataID fieldToClone_;

}; // class CloneBlockDataHandling






template< typename Field_T >
class FlattenedShallowCopyBlockDataHandling : public blockforest::AlwaysInitializeBlockDataHandling< typename Field_T::FlattenedField >
{
public:

   FlattenedShallowCopyBlockDataHandling( const ConstBlockDataID & fieldToClone ) :
      fieldToClone_( fieldToClone )
   {}

   typename Field_T::FlattenedField * initialize( IBlock * const block )
   {
      const Field_T * toClone = block->template getData< Field_T >( fieldToClone_ );
      return toClone->flattenedShallowCopy();
   }

private:

   ConstBlockDataID fieldToClone_;

}; // class FlattenedShallowCopyBlockDataHandling



} // namespace field
} // namespace walberla
