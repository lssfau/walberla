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
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Block.h"

#include "domain_decomposition/BlockDataHandling.h"



namespace walberla {
namespace blockforest {



template< typename T >
class BlockDataHandling : public domain_decomposition::BlockDataHandling<T>
{
public:
   virtual ~BlockDataHandling() {}

   /// must be thread-safe !
   virtual void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child ) = 0;
   /// must be thread-safe !
   virtual void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) = 0;
   
   /// must be thread-safe !
   virtual T * deserializeCoarseToFine( Block * const block ) = 0;
   /// must be thread-safe !
   virtual T * deserializeFineToCoarse( Block * const block ) = 0;
   
   /// must be thread-safe !
   virtual void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) = 0;
   /// must be thread-safe !
   virtual void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child ) = 0;
};



template< typename T >
class AlwaysInitializeBlockDataHandling : public BlockDataHandling<T>
{
public:
   ~AlwaysInitializeBlockDataHandling() {}

   void serialize( IBlock * const, const BlockDataID &, mpi::SendBuffer & ) {}
   void serializeCoarseToFine( Block * const, const BlockDataID &, mpi::SendBuffer &, const uint_t ) {}
   void serializeFineToCoarse( Block * const, const BlockDataID &, mpi::SendBuffer & ) {}

   T * deserialize( IBlock * const block ) { return this->initialize( block ); }
   T * deserializeCoarseToFine( Block * const block ) { return this->initialize( block ); }
   T * deserializeFineToCoarse( Block * const block ) { return this->initialize( block ); }

   void deserialize( IBlock * const, const BlockDataID &, mpi::RecvBuffer & ) {}
   void deserializeCoarseToFine( Block * const, const BlockDataID &, mpi::RecvBuffer & ) {}
   void deserializeFineToCoarse( Block * const, const BlockDataID &, mpi::RecvBuffer &, const uint_t ) {}
};



template< typename T >
class AlwaysCreateBlockDataHandling : public AlwaysInitializeBlockDataHandling<T>
{
public:
   ~AlwaysCreateBlockDataHandling() {}

   T * initialize( IBlock * const /*block*/ ) {return new T();}
};



namespace internal {



class BlockDataHandlingWrapper : public domain_decomposition::internal::BlockDataHandlingWrapper
{
public:
   typedef domain_decomposition::internal::BlockData BlockData;

   virtual ~BlockDataHandlingWrapper() {}
   
   virtual void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child ) = 0;
   virtual void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) = 0;
   
   virtual BlockData * deserializeCoarseToFine( Block * const block ) = 0;
   virtual BlockData * deserializeFineToCoarse( Block * const block ) = 0;
   
   virtual void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) = 0;
   virtual void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child ) = 0;
};



template< typename T >
class BlockDataHandlingHelper : public BlockDataHandlingWrapper
{
public:
   typedef domain_decomposition::internal::BlockData BlockData;

   BlockDataHandlingHelper( const shared_ptr< BlockDataHandling<T> > & dataHandling ) : dataHandling_( dataHandling ) {}
   ~BlockDataHandlingHelper() {}
   
   BlockData * initialize( IBlock * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      T * ptr = dataHandling_->initialize( block );
      return ptr ? new BlockData( ptr ) : NULL;
   }
   
   void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->serialize( block, id, buffer );
   }
   
   void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->serializeCoarseToFine( block, id, buffer, child );
   }
   
   void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->serializeFineToCoarse( block, id, buffer );
   }
   
   BlockData * deserialize( IBlock * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      T * ptr = dataHandling_->deserialize( block );
      return ptr ? new BlockData( ptr ) : NULL;
   }
   
   BlockData * deserializeCoarseToFine( Block * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      T * ptr = dataHandling_->deserializeCoarseToFine( block );
      return ptr ? new BlockData( ptr ) : NULL;
   }
   
   BlockData * deserializeFineToCoarse( Block * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      T * ptr = dataHandling_->deserializeFineToCoarse( block );
      return ptr ? new BlockData( ptr ) : NULL;
   }
   
   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->deserialize( block, id, buffer );
   }
   
   void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->deserializeCoarseToFine( block, id, buffer );
   }   
   
   void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      dataHandling_->deserializeFineToCoarse( block, id, buffer, child );
   }
   
private:

   shared_ptr< BlockDataHandling<T> > dataHandling_; 
};



} // namespace internal
} // namespace blockforest



namespace domain_decomposition {
class BlockStorage;
}

namespace blockforest {
namespace internal {

class BlockDataHandlingAdder
{
public:

   BlockDataHandlingAdder( domain_decomposition::BlockStorage & storage, const std::string & identifier = std::string() ) :
      storage_( storage ), identifier_( identifier )
   {}

   template< typename T >
   BlockDataHandlingAdder & operator<<( const domain_decomposition::BlockDataCreator<T> & bdc )
   {
      auto downcast = boost::dynamic_pointer_cast< blockforest::BlockDataHandling<T> >( bdc.dataHandling_ );
      if( downcast )
      {
         dataHandling_.add( make_shared< blockforest::internal::BlockDataHandlingHelper<T> >( downcast ), bdc.requiredSelectors_,
                            bdc.incompatibleSelectors_, bdc.identifier_ );
      }
      else
      {
         dataHandling_.add( make_shared< domain_decomposition::internal::BlockDataHandlingHelper<T> >( bdc.dataHandling_ ), bdc.requiredSelectors_,
                            bdc.incompatibleSelectors_, bdc.identifier_ );
      }

      return *this;
   }

   operator BlockDataID();

private:

   domain_decomposition::BlockStorage & storage_;
   std::string identifier_;
   domain_decomposition::internal::SelectableBlockDataHandlingWrapper dataHandling_;
};



} // namespace internal
} // namespace blockforest
} // namespace walberla
