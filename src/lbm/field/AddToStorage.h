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
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "PdfField.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/uid/SUID.h"
#include "field/blockforest/BlockDataHandling.h"

namespace walberla {
namespace lbm {



namespace internal {
   
template< typename LatticeModel_T >
class PdfFieldHandling : public field::BlockDataHandling< PdfField<LatticeModel_T>,
                                                          LatticeModel_T::Stencil::D == 2 >
{
public:

   typedef PdfField<LatticeModel_T> PdfField_T;
   typedef field::BlockDataHandling< PdfField_T, LatticeModel_T::Stencil::D == 2 > Base_T;

   PdfFieldHandling( const weak_ptr< StructuredBlockStorage > & blocks, const LatticeModel_T & latticeModel,
                     const bool _initialize, const Vector3<real_t> & initialVelocity, const real_t initialDensity,
                     const uint_t nrOfGhostLayers, const field::Layout & layout ) :
      blocks_( blocks ), latticeModel_( latticeModel ),
      initialize_( _initialize ), initialVelocity_( initialVelocity ), initialDensity_( initialDensity ),
      nrOfGhostLayers_( nrOfGhostLayers ), layout_( layout ) {}

   inline void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
   {
      packLatticeModel( block, id, buffer );
      Base_T::serialize( block, id, buffer );
   }

   void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child )
   {
      packLatticeModel( block, id, buffer );
      Base_T::serializeCoarseToFine( block, id, buffer, child );
   }

   void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
   {
      packLatticeModel( block, id, buffer );
      Base_T::serializeFineToCoarse( block, id, buffer );
   }

   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
   {
      unpackLatticeModel( block, id, buffer );
      Base_T::deserialize( block, id, buffer );
   }

   void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
   {
      unpackLatticeModel( block, id, buffer );
      Base_T::deserializeCoarseToFine( block, id, buffer );
   }

   void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child )
   {
      unpackLatticeModel( block, id, buffer );
      Base_T::deserializeFineToCoarse( block, id, buffer, child );
   }

protected:

   PdfField<LatticeModel_T> * allocate( IBlock * const block )
   {
      return allocateDispatch( block, initialize_, initialDensity_ );
   }

   PdfField<LatticeModel_T> * reallocate( IBlock * const block )
   {
#ifdef NDEBUG
      return allocateDispatch( block, false, initialDensity_ );
#else
      return allocateDispatch( block, true, std::numeric_limits<real_t>::quiet_NaN() );
#endif
   }

private:

   void packLatticeModel( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) const
   {
      const PdfField_T * field = block->template getData< PdfField_T >(id);
      WALBERLA_CHECK_NOT_NULLPTR( field );
      buffer << field->latticeModel();
   }

   void unpackLatticeModel( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) const
   {
      PdfField_T * field = block->template getData< PdfField_T >(id);
      WALBERLA_CHECK_NOT_NULLPTR( field );

      LatticeModel_T latticeModel = field->latticeModel();
      buffer >> latticeModel;

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks );

      latticeModel.configure( *block, *blocks );
      field->resetLatticeModel( latticeModel );
   }

   PdfField<LatticeModel_T> * allocateDispatch( IBlock * const block, const bool _initialize, const real_t initialDensity )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks );

      latticeModel_.configure( *block, *blocks );

      return new PdfField_T( blocks->getNumberOfXCells( *block ), blocks->getNumberOfYCells( *block ), blocks->getNumberOfZCells( *block ),
                             latticeModel_, _initialize, initialVelocity_, initialDensity, nrOfGhostLayers_, layout_ );
   }

   weak_ptr< StructuredBlockStorage > blocks_;

   LatticeModel_T    latticeModel_;
   bool              initialize_;
   Vector3< real_t > initialVelocity_;
   real_t            initialDensity_;
   uint_t            nrOfGhostLayers_;
   field::Layout     layout_;

}; // class PdfFieldHandling

} // namespace internal



template< typename LatticeModel_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                  const LatticeModel_T & latticeModel,
                                  const field::Layout & layout = field::zyxf,
                                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                   blocks, latticeModel, true, Vector3<real_t>(0), real_t(1), uint_t(1), layout ),
                                identifier, requiredSelectors, incompatibleSelectors );
}



template< typename LatticeModel_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                  const LatticeModel_T & latticeModel,
                                  const uint_t ghostLayers,
                                  const field::Layout & layout = field::zyxf,
                                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                   blocks, latticeModel, true, Vector3<real_t>(0), real_t(1), ghostLayers, layout ),
                                identifier, requiredSelectors, incompatibleSelectors );
}



template< typename LatticeModel_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                  const LatticeModel_T & latticeModel,
                                  const Vector3< real_t > & initialVelocity, const real_t initialDensity,
                                  const field::Layout & layout = field::zyxf,
                                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                   blocks, latticeModel, true, initialVelocity, initialDensity, uint_t(1), layout ),
                                identifier, requiredSelectors, incompatibleSelectors );
}



template< typename LatticeModel_T, typename BlockStorage_T >
BlockDataID addPdfFieldToStorage( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                                  const LatticeModel_T & latticeModel,
                                  const Vector3< real_t > & initialVelocity, const real_t initialDensity,
                                  const uint_t ghostLayers,
                                  const field::Layout & layout = field::zyxf,
                                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                   blocks, latticeModel, true, initialVelocity, initialDensity, ghostLayers, layout ),
                                identifier, requiredSelectors, incompatibleSelectors );
}






template< typename LatticeModel_T >
struct PdfFieldCreator : public domain_decomposition::BlockDataCreator< lbm::PdfField< LatticeModel_T > >
{
   PdfFieldCreator( const shared_ptr< StructuredBlockStorage > & blocks,
                    const std::string & identifier, const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                    const LatticeModel_T & latticeModel,
                    const field::Layout & layout = field::zyxf ) :
      domain_decomposition::BlockDataCreator< lbm::PdfField< LatticeModel_T > >( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                                                                    blocks, latticeModel, false, Vector3<real_t>(0), real_t(1), uint_t(1), layout ),
                                                                                 identifier, requiredSelectors, incompatibleSelectors )
   {}

   PdfFieldCreator( const shared_ptr< StructuredBlockStorage > & blocks,
                    const std::string & identifier, const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                    const LatticeModel_T & latticeModel, const uint_t ghostLayers,
                    const field::Layout & layout = field::zyxf ) :
      domain_decomposition::BlockDataCreator< lbm::PdfField< LatticeModel_T > >( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                                                                    blocks, latticeModel, false, Vector3<real_t>(0), real_t(1), ghostLayers, layout ),
                                                                                 identifier, requiredSelectors, incompatibleSelectors )
   {}

   PdfFieldCreator( const shared_ptr< StructuredBlockStorage > & blocks,
                    const std::string & identifier, const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                    const LatticeModel_T & latticeModel, const Vector3< real_t > & initialVelocity, const real_t initialDensity,
                    const field::Layout & layout = field::zyxf ) :
      domain_decomposition::BlockDataCreator< lbm::PdfField< LatticeModel_T > >( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                                                                    blocks, latticeModel, true, initialVelocity, initialDensity, uint_t(1), layout ),
                                                                                 identifier, requiredSelectors, incompatibleSelectors )
   {}

   PdfFieldCreator( const shared_ptr< StructuredBlockStorage > & blocks,
                    const std::string & identifier, const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                    const LatticeModel_T & latticeModel, const Vector3< real_t > & initialVelocity, const real_t initialDensity, const uint_t ghostLayers,
                    const field::Layout & layout = field::zyxf ) :
      domain_decomposition::BlockDataCreator< lbm::PdfField< LatticeModel_T > >( make_shared< internal::PdfFieldHandling< LatticeModel_T > >(
                                                                                    blocks, latticeModel, true, initialVelocity, initialDensity, ghostLayers, layout ),
                                                                                 identifier, requiredSelectors, incompatibleSelectors )
   {}
};



} // namespace lbm
} // namespace walberla
