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
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockDataID.h"
#include "IBlock.h"

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "core/selectable/SetSelectableObject.h"
#include "core/uid/GlobalState.h"
#include "core/uid/SUID.h"

#include <string>



namespace walberla {
namespace domain_decomposition {



template< typename T >
class BlockDataHandling
{
public:
   using value_type = T;

   virtual ~BlockDataHandling() = default;
   
   /// must be thread-safe !
   virtual T * initialize( IBlock * const block ) = 0;
   
   /// must be thread-safe !
   virtual void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) = 0;
   /// must be thread-safe !
   virtual T * deserialize( IBlock * const block ) = 0;
   /// must be thread-safe !
   virtual void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) = 0;
};



template< typename T >
class AlwaysInitializeBlockDataHandling : public BlockDataHandling<T>
{
public:
   ~AlwaysInitializeBlockDataHandling() override = default;

   void serialize( IBlock * const, const BlockDataID &, mpi::SendBuffer & ) override {}
   T * deserialize( IBlock * const block ) override { return this->initialize( block ); }
   void deserialize( IBlock * const, const BlockDataID &, mpi::RecvBuffer & ) override {}
};



namespace internal {

template< typename T >
class BlockDataHandlingFunctionAdaptor : public BlockDataHandling<T>
{
public:

   using Function = std::function<T *(IBlock *const)>;

   BlockDataHandlingFunctionAdaptor( const Function & function ) : function_( function ) {}

   ~BlockDataHandlingFunctionAdaptor() override = default;

   T * initialize( IBlock * const block ) override { return function_( block ); }
   
   void serialize( IBlock * const, const BlockDataID &, mpi::SendBuffer & ) override
   {
      WALBERLA_ABORT( "You are trying to serialize a block data item for which only an initialization function was registered" )
#ifdef __IBMCPP__
      return nullptr; // never reached, helps to suppress a warning from the IBM compiler
#endif
   }
   T * deserialize( IBlock * const ) override
   {
      WALBERLA_ABORT( "You are trying to deserialize a block data item for which only an initialization function was registered" )
#ifdef __IBMCPP__
      return nullptr; // never reached, helps to suppress a warning from the IBM compiler
#endif
   }
   void deserialize( IBlock * const, const BlockDataID &, mpi::RecvBuffer & ) override
   {
      WALBERLA_ABORT( "You are trying to deserialize a block data item for which only an initialization function was registered" )
   }

private:

   Function function_;
};

} // namespace internal



template< typename T >
struct BlockDataCreator
{
   template< typename U >
   BlockDataCreator( const shared_ptr< U > & dataHandling,
                     const std::string & identifier            = std::string(),
                     const Set<SUID> &   requiredSelectors     = Set<SUID>::emptySet(),
                     const Set<SUID> &   incompatibleSelectors = Set<SUID>::emptySet() ) :
      dataHandling_( dataHandling ), identifier_( identifier ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   BlockDataCreator( const std::function< T* ( IBlock * const block ) > & function,
                     const std::string & identifier            = std::string(),
                     const Set<SUID> &   requiredSelectors     = Set<SUID>::emptySet(),
                     const Set<SUID> &   incompatibleSelectors = Set<SUID>::emptySet() ) :
      dataHandling_( walberla::make_shared< internal::BlockDataHandlingFunctionAdaptor<T> >( function ) ), identifier_( identifier ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   shared_ptr< BlockDataHandling<T> > dataHandling_;

   std::string identifier_;
   Set<SUID>   requiredSelectors_;
   Set<SUID>   incompatibleSelectors_;
};



class BlockStorage;
namespace internal {



class BlockDataHandlingWrapper
{
public:
   virtual ~BlockDataHandlingWrapper() = default;

   virtual BlockData * initialize( IBlock * const block ) = 0;
   
   virtual void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) = 0;
   virtual BlockData * deserialize( IBlock * const block ) = 0;
   virtual void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) = 0;
};



template< typename T >
class BlockDataHandlingHelper : public BlockDataHandlingWrapper
{
public:

   BlockDataHandlingHelper( const shared_ptr< BlockDataHandling<T> > & dataHandling ) : dataHandling_( dataHandling ) {}
  
   BlockData * initialize( IBlock * const block ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block )
      T * ptr = dataHandling_->initialize( block );
      return ptr ? new BlockData( ptr ) : nullptr;
   }
   
   void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block )
      dataHandling_->serialize( block, id, buffer );
   }
   
   BlockData * deserialize( IBlock * const block ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block )
      T * ptr = dataHandling_->deserialize( block );
      return ptr ? new BlockData( ptr ) : nullptr;
   }
   
   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block )
      dataHandling_->deserialize( block, id, buffer );
   }   
   
private:

   shared_ptr< BlockDataHandling<T> > dataHandling_; 
};



using SelectableBlockDataHandlingWrapper = selectable::SetSelectableObject<shared_ptr<BlockDataHandlingWrapper>, SUID>;



class BlockDataHandlingAdder
{
public:

   BlockDataHandlingAdder( BlockStorage & storage, const std::string & identifier = std::string() ) :
      storage_( storage ), identifier_( identifier )
   {}

   template< typename T >
   BlockDataHandlingAdder & operator<<( const BlockDataCreator<T> & bdc )
   {
      dataHandling_.add( walberla::make_shared< BlockDataHandlingHelper<T> >( bdc.dataHandling_ ), bdc.requiredSelectors_,
                         bdc.incompatibleSelectors_, bdc.identifier_ );
      return *this;
   }

   operator BlockDataID();

private:

   BlockStorage & storage_;
   std::string identifier_;
   SelectableBlockDataHandlingWrapper dataHandling_;
};



class BlockDataItem
{
public:
   
   BlockDataItem( const BlockDataID & id, const std::string & identifier, const SelectableBlockDataHandlingWrapper & dataHandling ) :
      id_( id ), identifier_( identifier ), dataHandling_( dataHandling )
   {}
   
   bool operator==( const BlockDataItem & rhs ) const
   {
      return id_ == rhs.id_ && identifier_ == rhs.identifier_; // TODO: include comparison of 'dataHandling_'
   }
   bool operator!=( const BlockDataItem & rhs ) const { return !operator==( rhs ); }

   const BlockDataID & getId() const { return id_; }
   const std::string & getIdentifier() const { return identifier_; }
   
   shared_ptr< BlockDataHandlingWrapper > getDataHandling( IBlock const * const block, const Set<SUID> & state = Set<SUID>::emptySet() )
   {
      shared_ptr< BlockDataHandlingWrapper > dataHandling;
      
      Set<SUID> selection( uid::globalState() + block->getState() + state );
      size_t numMatches = dataHandling_.get( dataHandling, selection );

      if( numMatches > size_t(1) )
      {
         WALBERLA_ABORT( "More than one data handling object found for block data \"" << identifier_ << "\"\n"
                         " - number of matching data handling objects: " << numMatches << "\n"
                         " - block ID: " << block->getId() << "\n"
                         " - block state: " << block->getState() << "\n"
                         " - global state: " << uid::globalState() << "\n"
                         " - additional state: " << state << "\n" 
                         " - \"selector\": " << selection )
      }
      
      return dataHandling;
   }

private:

   BlockDataID id_;
   std::string identifier_;
   SelectableBlockDataHandlingWrapper dataHandling_;
};



} // namespace internal
} // namespace domain_decomposition
} // namespace walberla
