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
//! \file AddToStorage.h
//! \ingroup lbm_generated
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "PdfField.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/uid/SUID.h"
#include "field/blockforest/BlockDataHandling.h"

namespace walberla::lbm_generated {

namespace internal {
   
template< typename LatticeStorageSpecification_T >
class PdfFieldHandling : public field::BlockDataHandling< PdfField<LatticeStorageSpecification_T>,
                                                          LatticeStorageSpecification_T::Stencil::D == 2 >
{
public:

   using PdfField_T = PdfField<LatticeStorageSpecification_T>;
   using Base_T = field::BlockDataHandling<PdfField_T, LatticeStorageSpecification_T::Stencil::D == 2>;

   PdfFieldHandling( const weak_ptr< StructuredBlockStorage > & blocks, const LatticeStorageSpecification_T & storageSpecification,
                     const uint_t nrOfGhostLayers, const field::Layout & layout, const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr ) :
      blocks_( blocks ), storageSpecification_( storageSpecification ),
      nrOfGhostLayers_( nrOfGhostLayers ), layout_( layout ), alloc_( alloc ){}

   inline void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override
   {
      Base_T::serialize( block, id, buffer );
   }

   void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child ) override
   {
      Base_T::serializeCoarseToFine( block, id, buffer, child );
   }

   void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override
   {
      Base_T::serializeFineToCoarse( block, id, buffer );
   }

   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override
   {
      Base_T::deserialize( block, id, buffer );
   }

   void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override
   {
      Base_T::deserializeCoarseToFine( block, id, buffer );
   }

   void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child ) override
   {
      Base_T::deserializeFineToCoarse( block, id, buffer, child );
   }

protected:

   PdfField<LatticeStorageSpecification_T> * allocate( IBlock * const block ) override
   {
      return allocateDispatch( block );
   }

   PdfField<LatticeStorageSpecification_T> * reallocate( IBlock * const block ) override
   {
      return allocateDispatch( block );
   }

private:


   PdfField<LatticeStorageSpecification_T> * allocateDispatch( IBlock * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block )

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks )

      return new PdfField_T( blocks->getNumberOfXCells( *block ), blocks->getNumberOfYCells( *block ), blocks->getNumberOfZCells( *block ),
                            storageSpecification_, nrOfGhostLayers_, layout_, alloc_ );
   }

   weak_ptr< StructuredBlockStorage > blocks_;
   LatticeStorageSpecification_T    storageSpecification_;

   uint_t            nrOfGhostLayers_;
   field::Layout     layout_;
   shared_ptr< field::FieldAllocator<real_t> > alloc_;

}; // class PdfFieldHandling

} // namespace internal



template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                  const LatticeStorageSpecification_T & storageSpecification,
                                  const uint_t ghostLayers,
                                  const field::Layout & layout = field::fzyx,
                                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
                                  const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr)
{
   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                   blocks, storageSpecification, ghostLayers, layout, alloc ),
                                identifier, requiredSelectors, incompatibleSelectors );
}

template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                 const LatticeStorageSpecification_T & storageSpecification,
                                 const field::Layout & layout = field::fzyx,
                                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
                                 const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr)
{
   auto ghostLayers = uint_c(1);

   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                  blocks, storageSpecification, ghostLayers, layout, alloc ),
                               identifier, requiredSelectors, incompatibleSelectors );
}

template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                 const LatticeStorageSpecification_T & storageSpecification,
                                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
                                 const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr)
{
   auto ghostLayers = uint_c(1);
   auto layout = field::fzyx;

   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                  blocks, storageSpecification, ghostLayers, layout, alloc ),
                               identifier, requiredSelectors, incompatibleSelectors );
}

template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                 const LatticeStorageSpecification_T & storageSpecification,
                                 const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr)
{
   auto ghostLayers = uint_c(1);
   auto layout = field::fzyx;
   auto requiredSelectors = Set<SUID>::emptySet();
   auto incompatibleSelectors = Set<SUID>::emptySet();

   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                  blocks, storageSpecification, ghostLayers, layout, alloc ),
                               identifier, requiredSelectors, incompatibleSelectors );
}

template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                 const LatticeStorageSpecification_T & storageSpecification,
                                 const field::Layout & layout = field::fzyx,
                                 const shared_ptr< field::FieldAllocator<real_t> > alloc = nullptr)
{
   auto ghostLayers = uint_c(1);
   auto requiredSelectors = Set<SUID>::emptySet();
   auto incompatibleSelectors = Set<SUID>::emptySet();

   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                  blocks, storageSpecification, ghostLayers, layout, alloc ),
                               identifier, requiredSelectors, incompatibleSelectors );
}

template< typename LatticeStorageSpecification_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                 const LatticeStorageSpecification_T & storageSpecification,
                                 const uint_t ghostLayers,
                                 const field::Layout & layout,
                                 const shared_ptr< field::FieldAllocator<real_t> > alloc)
{
   auto requiredSelectors = Set<SUID>::emptySet();
   auto incompatibleSelectors = Set<SUID>::emptySet();

   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeStorageSpecification_T > >(
                                  blocks, storageSpecification, ghostLayers, layout, alloc ),
                               identifier, requiredSelectors, incompatibleSelectors );
}


} // namespace walberla::lbm_generated
