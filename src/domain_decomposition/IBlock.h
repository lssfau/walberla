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
//! \file IBlock.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockDataID.h"
#include "IBlockID.h"

#include "core/Abort.h"
#include "core/debug/demangle.h"
#include "core/NonCopyable.h"
#include "core/debug/Debug.h"
#include "core/math/AABB.h"
#include "core/uid/SUID.h"

#include <typeinfo>
#include <vector>


namespace walberla {
namespace domain_decomposition {



namespace internal {

/// wrapper class for any kind of block data (used only internally, never to be seen in the interface of public member functions)
/// see: http://www.drdobbs.com/cpp/twisting-the-rtti-system-for-safe-dynami/229401004#
class BlockData : private NonCopyable
{
private:

   class DataBase {
   public:
      virtual ~DataBase() = default;
      virtual bool operator==( const DataBase& rhs ) const = 0;
      virtual bool isOfSameType( const DataBase& rhs ) const = 0;
   };

   template< typename T >
   class Data : public DataBase {
   public:
      Data( T* data ) : data_( data ) {}
      ~Data() override { delete data_; }
      bool operator==( const DataBase& rhs ) const override {
         const Data<T>* rhsData = dynamic_cast< const Data<T>* >( &rhs );
         return ( rhsData == &rhs ) && ( *data_ == *(rhsData->data_) ); // every object that is registered as block data
                                                                        // must be comparable with "==" !
      }
      bool isOfSameType( const DataBase& rhs ) const override { return dynamic_cast< const Data<T>* >( &rhs ) == &rhs; }
   private:
      T* data_;
   };

public:

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning( disable : 4670 )
#  pragma warning( disable : 4673 )
#endif //_MSC_VER

   template< typename T >
   BlockData( T* ptr ) : ptr_( ptr ), thr_( &thrower<T> )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( ptr );
      data_ = new Data<T>( ptr );
#ifndef NDEBUG
      typeInfo_ = typeid(T).name();
#endif
   }

#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER

   ~BlockData() { delete data_; }

   bool operator==( const BlockData& rhs ) const { return *data_ == *rhs.data_; } // every object that is registered as block data
   bool operator!=( const BlockData& rhs ) const { return !operator==( rhs ); }   // must be comparable with "==" !

   template< typename U >
   const U* get() const {
      try { thr_(ptr_); }
      catch ( U* ptr ) { return ptr; }
      catch (...) {}
#ifndef NDEBUG
      WALBERLA_ABORT( "BlockData access type violation! (The block data you added is of a different type than the block data you are trying to access!)"
                      "\nThe original data type was:       " << debug::demangle( typeInfo_ ) <<
                      "\nYou try to retrieve data of type: " << debug::demangle( typeid(U).name() ) )
#else
      WALBERLA_ABORT( "BlockData access type violation! (The block data you added is of a different type than the block data you are trying to access!)" )
#endif
#ifdef __IBMCPP__
      return nullptr; // never reached, helps to suppress a warning from the IBM compiler
#endif
   }

   template< typename U >
   U* get() { return const_cast< U* >( static_cast< const BlockData* >(this)->template get<U>() ); }

   /// Attention: Calling this function ONLY works if type U is identical to the type T that was once added.
   ///            Retrieving a pointer to a base class of class T that was once added does not work!
   template< typename U >
   const U* uncheckedFastGet() const { return static_cast< const U* >( ptr_ ); }

   /// Attention: Calling this function ONLY works if type U is identical to the type T that was once added.
   ///            Retrieving a pointer to a base class of class T that was once added does not work!
   template< typename U >
         U* uncheckedFastGet()       { return const_cast< U* >( static_cast< const BlockData* >(this)->template uncheckedFastGet<U>() ); }

   template< typename U >
   bool isOfType() const {
      return dynamic_cast< const Data<U>* >( data_ ) == data_;
   }

   template< typename U >
   bool isClassOrSubclassOf() const {
      try { thr_(ptr_); }
      catch ( U* ) { return true; }
      catch (...) {}
      return false;
   }

   template< typename U >
   bool isSubclassOf() const {
      return isClassOrSubclassOf<U>() && !isOfType<U>();
   }

   bool isOfSameType( const BlockData& bd ) const {
      return data_->isOfSameType( *bd.data_ );
   }

private:

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning( disable : 4670 )
#  pragma warning( disable : 4673 )
#endif //_MSC_VER
   template< typename T > static void thrower( void* ptr ) { throw static_cast< T* >( ptr ); }
#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER

   DataBase* data_;
   void* ptr_;
   void (*thr_)( void* );

#ifndef NDEBUG
   std::string typeInfo_; // used only in debug mode!
#endif
};

} // namespace internal



class           BlockStorage; // forward declaration
class StructuredBlockStorage; // forward declaration



//**********************************************************************************************************************
/*!
*   \brief Base class for blocks (blocks are used to partition the simulation space: blocks are rectangular parts of the
*          simulation space that manage all the data that is assigned to their part of the simulation space)
*
*   The class 'IBlock' already provides some basic functionality that every type of block will possess - regardless of
*   the actual implementation of any derived block class. Every block ...
*
*      - ... has a well defined state - that is a set of SUIDs that describes certain attributes related to this block.
*        The state of a block can be used to determine what kind of data is assigned to this block and what kind of
*        operations/computations are performed with this block.
*      - ... possesses an axis-aligned bounding box that defines the exact part of the simulation space that is assigned
*        to this block.
*      - ... can hold any kind of data, which means every block can hold any number of dynamically allocated objects of
*        any type. This "block data" is assigned via block data initialization functions which are managed by the block
*        storage data structure that governs this block. For more information on this procedure refer to 'BlockStorage'.
*        Data stored within a block can be retrieved by calling the member functions 'getData'.
*
*   'IBlock' also acts as an interface class which requires every derived class to implement a member function 'getId()'
*   that returns a block ID of type 'IBlockID' (see classes 'IBlockID' and 'BlockStorage') which can be used to uniquely
*   identify a block within its governing block storage data structure.
*/
//**********************************************************************************************************************

class IBlock : private NonCopyable
{
public:

   using BlockData = internal::BlockData;

   friend class           BlockStorage;
   friend class StructuredBlockStorage;

   virtual const IBlockID& getId() const = 0;

   bool operator==( const IBlock& rhs ) const;
   bool operator!=( const IBlock& rhs ) const { return !operator==( rhs ); }

   const Set<SUID>& getState() const { return state_; }
   void setState( const Set<SUID>& state ) { state_  = state; } // Changing the state of a block only changes the internal 'state_' variable.
   void addState( const Set<SUID>& state ) { state_ += state; } // If other blocks on the same or on another process store information about the
   void addState( const SUID&      state ) { state_ += state; } // state of this block, this information gets invalidated! In order to guarantee
   void clearState() { state_.clear(); }                        // consistency, block state changes must somehow be communicated.

   const AABB& getAABB() const { return aabb_; }

   inline bool isBlockDataAllocated( const ConstBlockDataID & index ) const;

   template< typename T >
   inline const T* getData( const ConstBlockDataID & index ) const;

   template< typename T >
   inline const T* getData( const BlockDataID & index ) const;
   template< typename T >
   inline       T* getData( const BlockDataID & index );

   inline void deleteData( const BlockDataID & index );

   template< typename T >
   inline bool isDataOfType( const ConstBlockDataID & index ) const;

   template< typename T >
   inline bool isDataClassOrSubclassOf( const ConstBlockDataID & index ) const;

   template< typename T >
   inline bool isDataSubclassOf( const ConstBlockDataID & index ) const;

   inline bool isDataOfSameType( const ConstBlockDataID & indexA, const ConstBlockDataID & indexB ) const;

   const BlockStorage& getBlockStorage() const { return storage_; } ///< returns a reference to the governing block storage data structure
         BlockStorage& getBlockStorage()       { return storage_; } ///< returns a reference to the governing block storage data structure

   template< typename T >
   inline const T* uncheckedFastGetData( const ConstBlockDataID & index ) const;
   template< typename T >
   inline const T* uncheckedFastGetData( const BlockDataID & index ) const;
   template< typename T >
   inline       T* uncheckedFastGetData( const BlockDataID & index );
         
protected:

   IBlock( BlockStorage& storage, const AABB& aabb, const IBlockID::IDType& id, const bool local = true );

   virtual ~IBlock(); ///< Must not be made public! No one should be allowed to delete a variable of type 'IBlock*'

   virtual bool equal( const IBlock* rhs ) const = 0;



   AABB aabb_; ///< an axis-aligned bounding box that defines the exact part of the simulation space that is assigned to this block

private:

   IBlock(); ///< Must not be made public or protected! Derived classes must call one of the available public/protected constructors.

   /// helper function for assigning data to this block (called by 'BlockStorage'), must not be made public or protected!
   inline void addData( const BlockDataID & index, BlockData* const data );

   Set<SUID> state_; ///< a set of SUIDs the describes the current state of this block

   std::vector< BlockData* > data_; ///< the data assigned to this block

   BlockStorage& storage_; ///< governing block storage data structure
   int           storageIndex_; ///< used internally by 'IBlock' and 'BlockStorage' in order to (quite easily) allow iterators for 'BlockStorage'

}; // class IBlock



//**********************************************************************************************************************
/*!
*   Returns true if block data is allocated on this block at index "index"
*/
//**********************************************************************************************************************
inline bool IBlock::isBlockDataAllocated( const ConstBlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   return ( data_[index] != nullptr );
}



//**********************************************************************************************************************
/*!
*   For documentation of this function see "const T* IBlock::getData( const BlockDataID index ) const".
*/
//**********************************************************************************************************************
template< typename T >
inline const T* IBlock::getData( const ConstBlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == nullptr )
      return nullptr;

   return data_[index]->template get< T >();
}



//**********************************************************************************************************************
/*!
*   Function for retrieving data - data that was previously added by the governing block storage data structure class
*   (see class 'BlockStorage') via block data initialization functions.
*   Data cannot be added manually by for example calling 'addData', but must be added through the governing block
*   storage data structure that keeps track of assigning BlockDataIDs (BlockDataIDs are required for
*   accessing/retrieving data stored within blocks).
*   Retrieving data will fail if T is neither the same type nor a base class of the data that was added.
*/
//**********************************************************************************************************************
template< typename T >
inline const T* IBlock::getData( const BlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == nullptr )
      return nullptr;

   return data_[index]->template get< T >();
}



//**********************************************************************************************************************
/*!
*   For documentation of this function see "const T* IBlock::getData( const BlockDataID index ) const".
*/
//**********************************************************************************************************************
template< typename T >
inline T* IBlock::getData( const BlockDataID & index ) {

   return const_cast< T* >( static_cast< const IBlock* >( this )->getData<T>( index ) );
}



//**********************************************************************************************************************
/*!
*   Function for removing all data that corresponds to block data ID 'index'.
*   Further calls to "getData" with 'index' will return NULL.
*/
//**********************************************************************************************************************
inline void IBlock::deleteData( const BlockDataID & index )
{
   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );
   delete data_[index];
   data_[index] = nullptr;
}



//**********************************************************************************************************************
/*!
*   Returns true if and only if the data stored at position "index" is of type T. This function will return false if the
*   type of the data is a subclass of T! For polymorphic type checking see member functions "isDataClassOrSubclassOf"
*   and "isDataSubclassOf".
*/
//**********************************************************************************************************************
template< typename T >
inline bool IBlock::isDataOfType( const ConstBlockDataID & index ) const
{
   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );
   WALBERLA_ASSERT_NOT_NULLPTR( data_[index] );

   return data_[index]->template isOfType<T>();
}



//**********************************************************************************************************************
/*!
*   Returns true if the data stored at position "index" is of type T or if it is a subclass of T.
*/
//**********************************************************************************************************************
template< typename T >
inline bool IBlock::isDataClassOrSubclassOf( const ConstBlockDataID & index ) const
{
   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );
   WALBERLA_ASSERT_NOT_NULLPTR( data_[index] );

   return data_[index]->template isClassOrSubclassOf<T>();
}



//**********************************************************************************************************************
/*!
*   Returns true if and only if the type of the data stored at position "index" is a subclass of T. This function will
*   return false if the type of the data stored at position "index" is equal to T!
*/
//**********************************************************************************************************************
template< typename T >
inline bool IBlock::isDataSubclassOf( const ConstBlockDataID & index ) const
{
   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );
   WALBERLA_ASSERT_NOT_NULLPTR( data_[index] );

   return data_[index]->template isSubclassOf<T>();
}



//**********************************************************************************************************************
/*!
*   Returns true if the type of the data stored at position "indexA" is equal to the type of the data stored at
*   position "indexB".
*/
//**********************************************************************************************************************
inline bool IBlock::isDataOfSameType( const ConstBlockDataID & indexA, const ConstBlockDataID & indexB ) const
{
   WALBERLA_ASSERT_LESS( uint_t( indexA ), data_.size() );
   WALBERLA_ASSERT_LESS( uint_t( indexB ), data_.size() );
   WALBERLA_ASSERT_NOT_NULLPTR( data_[indexA] );
   WALBERLA_ASSERT_NOT_NULLPTR( data_[indexB] );

   return data_[indexA]->isOfSameType( *data_[indexB] );
}



inline void IBlock::addData( const BlockDataID & index, BlockData * const data )
{
   if( data_.size() <= index )
      data_.resize( index+1, nullptr );

   if( data != nullptr ) {
      WALBERLA_ASSERT_NULLPTR( data_[index] );
      data_[index] = data;
   }
}



//**********************************************************************************************************************
/*!
*   For documentation of this function see "const T* IBlock::uncheckedFastGetData( const BlockDataID index ) const".
*/
//**********************************************************************************************************************
template< typename T >
inline const T* IBlock::uncheckedFastGetData( const ConstBlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == nullptr )
      return nullptr;

   return data_[index]->template uncheckedFastGet< T >();
}



//**********************************************************************************************************************
/*!
*   Identical to "getData", however only works if type T is identical (!) to the type that was once added. If type T is
*   a base class of the class that was once added, this function will lead to undefined behavior of your program!
*   Use with extreme caution!
*/
//**********************************************************************************************************************
template< typename T >
inline const T* IBlock::uncheckedFastGetData( const BlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == nullptr )
      return nullptr;

   return data_[index]->template uncheckedFastGet< T >();
}



//**********************************************************************************************************************
/*!
*   For documentation of this function see "const T* IBlock::uncheckedFastGetData( const BlockDataID index ) const".
*/
//**********************************************************************************************************************
template< typename T >
inline T* IBlock::uncheckedFastGetData( const BlockDataID & index ) {

   return const_cast< T* >( static_cast< const IBlock* >( this )->uncheckedFastGetData<T>( index ) );
}



} // namespace domain_decomposition

using domain_decomposition::IBlock;

} // namespace walberla
