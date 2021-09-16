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
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Functions to add fields to blocks
//!
//! In waLBerla, a central data structure is the BlockStorage. A BlockStorage itself contains
//! blocks, which hold the actual data. To add data to a block one has to specify a so called 'BlockDataHandling'
//! that, among other things, is responsible for the allocation and initialization of the block data on each separate block.
//! To store, for example, a vector<int> at each block we have to write
//! a 'BlockDataHandling' that allocates a new vector<int> and returns a pointer to it. This 'BlockDataHandling'
//! is registered at the BlockStorage and called for each block. For details see documentation of
//! BlockStorage::addBlockData() and StructuredBlockStorage::addBlockData().
//! Often one wants to store fields as block data, so this file provides helper functions for this common task.
//! When we create fields as block data, we effectively introduce a new subdivision of the domain. The first
//! subdivision is the block structure, whereas the block is then again structured in cells.
//! The subdivision into cells is managed by a StructuredBlockStorage. The StructuredBlockStorage is a BlockStorage that
//! knows about the subdivision into cells like number of cells per block, mapping for local to global cells, etc.
//! The following functions can be used to add a FlagField or any GhostLayerField as a block data
//! to a StructuredBlockStorage. Additionally a BlockDataCreator for each field is provided, if one wants
//! to make the creation dependent on selectors
//
//======================================================================================================================

#pragma once

#include "FlagField.h"
#include "GhostLayerField.h"
#include "field/blockforest/BlockDataHandling.h"

#include <type_traits>


namespace walberla {
namespace field {



////////////////
// FLAG FIELD //
////////////////

template< typename FlagField_T, typename BlockStorage_T >
BlockDataID addFlagFieldToStorage( const shared_ptr< BlockStorage_T > & blocks,
                                   const std::string & identifier,
                                   const uint_t nrOfGhostLayers = uint_t(1),
                                   const bool alwaysInitialize = false,
                                   const std::function< void ( FlagField_T * field, IBlock * const block ) > & initFunction =
                                      std::function< void ( FlagField_T * field, IBlock * const block ) >(),
                                   const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   if( alwaysInitialize )
   {
      auto dataHandling = make_shared< field::AlwaysInitializeBlockDataHandling< FlagField_T > >( blocks, nrOfGhostLayers );
      dataHandling->addInitializationFunction( initFunction );
      return blocks->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors );
   }

   auto dataHandling = make_shared< field::DefaultBlockDataHandling< FlagField_T > >( blocks, nrOfGhostLayers );
   dataHandling->addInitializationFunction( initFunction );
   return blocks->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors );
}



template< typename FlagField_T, typename BlockStorage_T >
BlockDataID addFlagFieldToStorage( const shared_ptr< BlockStorage_T > & blocks,
                                   const std::string & identifier,
                                   const uint_t nrOfGhostLayers,
                                   const bool alwaysInitialize,
                                   const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return addFlagFieldToStorage< FlagField_T >( blocks, identifier, nrOfGhostLayers, alwaysInitialize,
                                                std::function< void ( FlagField_T * field, IBlock * const block ) >(),
                                                requiredSelectors, incompatibleSelectors );
}



///////////////////////
// GHOST LAYER FIELD //
///////////////////////

namespace internal {

template< typename GhostLayerField_T, typename BlockStorage_T, class Enable = void >
struct AddToStorage
{
   using Value_T = typename GhostLayerField_T::value_type;
   static BlockDataID add( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                           const typename GhostLayerField_T::value_type & initValue, const Layout layout, const uint_t nrOfGhostLayers,
                           const bool /*alwaysInitialize*/, const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction,
                           const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                           const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = defaultSize,
                           const shared_ptr< field::FieldAllocator<Value_T> > alloc = nullptr)
   {
      auto dataHandling = walberla::make_shared< field::AlwaysInitializeBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout, calculateSize, alloc );
      dataHandling->addInitializationFunction( initFunction );
      return blocks->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors );
   }
};

template< typename GhostLayerField_T, typename BlockStorage_T >
struct AddToStorage< GhostLayerField_T, BlockStorage_T,
                     typename std::enable_if< ( std::is_integral< typename GhostLayerField_T::value_type >::value || std::is_floating_point< typename GhostLayerField_T::value_type >::value ) &&
	                                            ! std::is_same< GhostLayerField_T, FlagField< typename GhostLayerField_T::value_type > >::value >::type >
{
   using Value_T = typename GhostLayerField_T::value_type;
   static BlockDataID add( const shared_ptr< BlockStorage_T > & blocks, const std::string & identifier,
                           const typename GhostLayerField_T::value_type & initValue, const Layout layout, const uint_t nrOfGhostLayers,
                           const bool alwaysInitialize, const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction,
                           const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors,
                           const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize = defaultSize,
                           const shared_ptr< field::FieldAllocator<Value_T> > alloc = nullptr)
   {
      if( alwaysInitialize )
      {
         auto dataHandling = walberla::make_shared< field::AlwaysInitializeBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout, calculateSize, alloc );
         dataHandling->addInitializationFunction( initFunction );
         return blocks->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors );
      }

      auto dataHandling = walberla::make_shared< field::DefaultBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout, calculateSize, alloc );
      dataHandling->addInitializationFunction( initFunction );
      return blocks->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors );
   }
};

} // namespace internal



template< typename GhostLayerField_T, typename BlockStorage_T >
BlockDataID addToStorage( const shared_ptr< BlockStorage_T > & blocks,
                          const std::string & identifier,
                          const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
                          const Layout layout = zyxf,
                          const uint_t nrOfGhostLayers = uint_t(1),
                          const bool alwaysInitialize = false,
                          const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction =
                             std::function< void ( GhostLayerField_T * field, IBlock * const block ) >(),
                          const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                          const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet())
{
   return internal::AddToStorage< GhostLayerField_T, BlockStorage_T >::add( blocks, identifier, initValue, layout, nrOfGhostLayers,
                                                                            alwaysInitialize, initFunction, requiredSelectors, incompatibleSelectors );
}



template< typename GhostLayerField_T, typename BlockStorage_T >
BlockDataID addToStorage( const shared_ptr< BlockStorage_T > & blocks,
                          const std::string & identifier,
                          const typename GhostLayerField_T::value_type & initValue,
                          const Layout layout,
                          const uint_t nrOfGhostLayers,
                          const bool alwaysInitialize,
                          const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return addToStorage< GhostLayerField_T >( blocks, identifier, initValue, layout, nrOfGhostLayers, alwaysInitialize,
                                             std::function< void ( GhostLayerField_T * field, IBlock * const block ) >(),
                                             requiredSelectors, incompatibleSelectors );
}



template< typename GhostLayerField_T, typename BlockStorage_T >
BlockDataID addToStorage( const shared_ptr< BlockStorage_T > & blocks,
                          const std::string & identifier,
                          const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) >& calculateSize,
                          const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
                          const Layout layout = zyxf,
                          const uint_t nrOfGhostLayers = uint_t(1),
                          const bool alwaysInitialize = false,
                          const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction =
                          std::function< void ( GhostLayerField_T * field, IBlock * const block ) >(),
                          const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                          const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
                          const shared_ptr< field::FieldAllocator<typename GhostLayerField_T::value_type> > alloc = nullptr)
{
   return internal::AddToStorage< GhostLayerField_T, BlockStorage_T >::add( blocks, identifier, initValue, layout, nrOfGhostLayers,
                                                                            alwaysInitialize, initFunction, requiredSelectors,
                                                                            incompatibleSelectors, calculateSize, alloc );
}


template< typename GhostLayerField_T, typename BlockStorage_T >
BlockDataID addToStorage( const shared_ptr< BlockStorage_T > & blocks,
                          const std::string & identifier,
                          const std::function< Vector3< uint_t > ( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) > calculateSize,
                          const typename GhostLayerField_T::value_type & initValue,
                          const Layout layout,
                          const uint_t nrOfGhostLayers,
                          const bool alwaysInitialize,
                          const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return addToStorage< GhostLayerField_T >( blocks, identifier, initValue, layout, nrOfGhostLayers, alwaysInitialize,
                                             std::function< void ( GhostLayerField_T * field, IBlock * const block ) >(),
                                             requiredSelectors, incompatibleSelectors, calculateSize );
}






//**********************************************************************************************************************
/*! Adds a copy of an existing field to BlockStorage
*
*   \tparam Field_T         the type of the field that should be cloned ( and the type that is created )
*   \tparam BlockStorage_T  the type of the BlockStorage ( will be deduced automatically )
*
*   \param blocks        BlockStorage where the original field is stored and the new one is created
*   \param fieldToClone  BlockDataID of the Field that is cloned
*   \param identifier    name for new the field ( displayed in GUI and debugging functions )
*/
//**********************************************************************************************************************
template< typename Field_T, typename BlockStorage_T >
BlockDataID addCloneToStorage( const shared_ptr< BlockStorage_T > & blocks,
                               ConstBlockDataID fieldToClone,
                               const std::string & identifier,
                               const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                               const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< field::CloneBlockDataHandling< Field_T > >( fieldToClone ),
                                identifier, requiredSelectors, incompatibleSelectors );
}







//**********************************************************************************************************************
/*! Adds a flattened shallow copy of an existing field to BlockStorage
*
*   \tparam Field_T         the type of the field that should be cloned and flattened
*   \tparam BlockStorage_T  the type of the BlockStorage ( will be deduced automatically )
*
*   \param blocks        BlockStorage where the original field is stored and the new one is created
*   \param fieldToClone  BlockDataID of the Field that is cloned
*   \param identifier    name for new the field ( displayed in GUI and debugging functions )
*/
//**********************************************************************************************************************
template< typename Field_T, typename BlockStorage_T >
BlockDataID addFlattenedShallowCopyToStorage( const shared_ptr< BlockStorage_T > & blocks,
                                              ConstBlockDataID fieldToClone,
                                              const std::string & identifier,
                                              const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                              const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< field::FlattenedShallowCopyBlockDataHandling< Field_T > >( fieldToClone ),
                                identifier, requiredSelectors, incompatibleSelectors );
}







//**********************************************************************************************************************
/*! BlockDataCreator for fields
*
* Use this class when different fields should be added for different selectors.
*
* Example for making the field layout selector dependent:
* \code
   BlockDataID fieldId = blocks->addBlockData( "my data field" )
      << field::Creator< MyField_T >( blocks, "my data field (zyxf)", NONE, FZYX, 0, field::zyxf )
      << field::Creator< MyField_T >( blocks, "my data field (fzyx)", FZYX, NONE, 0, field::fzyx );
* \endcode
*/
//**********************************************************************************************************************
template< typename GhostLayerField_T, class Enable = void >
struct Creator : public domain_decomposition::BlockDataCreator< GhostLayerField_T >
{
   Creator( const shared_ptr< StructuredBlockStorage > & blocks,
            const std::string & identifier,
            const Set<SUID> & requiredSelectors,
            const Set<SUID> & incompatibleSelectors,
            const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
            const Layout layout = zyxf,
            const uint_t nrOfGhostLayers = uint_t(1),
            const bool /*alwaysInitialize*/ = false,
            const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction =
               std::function< void ( GhostLayerField_T * field, IBlock * const block ) >() ) :
      domain_decomposition::BlockDataCreator< GhostLayerField_T >( shared_ptr< field::DefaultBlockDataHandling< GhostLayerField_T > >(),
                                                                   identifier, requiredSelectors, incompatibleSelectors )
   {
      auto dataHandling = make_shared< field::AlwaysInitializeBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout );
      dataHandling->addInitializationFunction( initFunction );
      this->dataHandling_ = dataHandling;

   }

   Creator( const std::string & identifier = std::string(),
            const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
            const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
            const Layout layout = zyxf,
            const uint_t nrOfGhostLayers = uint_t(1) ) :
      domain_decomposition::BlockDataCreator< GhostLayerField_T >( shared_ptr< DefaultBlockDataHandling< GhostLayerField_T > >(),
                                                                   identifier, requiredSelectors, incompatibleSelectors )
   {
      static_assert( never_true<GhostLayerField_T>::value,
                     "This constructor of class 'field::Creator' is not supported anymore! Please use the other, new constructor." );
   }
};

template< typename GhostLayerField_T >
struct Creator< GhostLayerField_T,
                typename std::enable_if< std::is_integral< typename GhostLayerField_T::value_type >::value ||
                                         std::is_floating_point< typename GhostLayerField_T::value_type >::value  >::type >
   : public domain_decomposition::BlockDataCreator< GhostLayerField_T >
{
   Creator( const shared_ptr< StructuredBlockStorage > & blocks,
            const std::string & identifier,
            const Set<SUID> & requiredSelectors,
            const Set<SUID> & incompatibleSelectors,
            const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
            const Layout layout = zyxf,
            const uint_t nrOfGhostLayers = uint_t(1),
            const bool alwaysInitialize = false,
            const std::function< void ( GhostLayerField_T * field, IBlock * const block ) > & initFunction =
               std::function< void ( GhostLayerField_T * field, IBlock * const block ) >() ) :
      domain_decomposition::BlockDataCreator< GhostLayerField_T >( shared_ptr< field::DefaultBlockDataHandling< GhostLayerField_T > >(),
                                                                   identifier, requiredSelectors, incompatibleSelectors )
   {
      if( alwaysInitialize )
      {
         auto dataHandling = make_shared< field::AlwaysInitializeBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout );
         dataHandling->addInitializationFunction( initFunction );
         this->dataHandling_ = dataHandling;
      }
      else
      {
         auto dataHandling = make_shared< field::DefaultBlockDataHandling< GhostLayerField_T > >( blocks, nrOfGhostLayers, initValue, layout );
         dataHandling->addInitializationFunction( initFunction );
         this->dataHandling_ = dataHandling;
      }
   }
   
   Creator( const std::string & identifier = std::string(),
            const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
            const typename GhostLayerField_T::value_type & initValue = typename GhostLayerField_T::value_type(),
            const Layout layout = zyxf,
            const uint_t nrOfGhostLayers = uint_t(1) ) :
      domain_decomposition::BlockDataCreator< GhostLayerField_T >( shared_ptr< DefaultBlockDataHandling< GhostLayerField_T > >(),
                                                                   identifier, requiredSelectors, incompatibleSelectors )
   {
      static_assert( never_true<GhostLayerField_T>::value,
                     "This constructor of class 'field::Creator' is not supported anymore! Please use the other, new constructor." );
   }
};



// specialization for class FlagField
template< typename T>
struct Creator< FlagField<T> > : public domain_decomposition::BlockDataCreator< FlagField<T> >
{
   Creator( const shared_ptr< StructuredBlockStorage > & blocks,
            const std::string & identifier,
            const Set<SUID> & requiredSelectors,
            const Set<SUID> & incompatibleSelectors,
            const uint_t nrOfGhostLayers = uint_t(1),
            const bool alwaysInitialize = false,
            const std::function< void ( FlagField<T> * field, IBlock * const block ) > & initFunction =
               std::function< void ( FlagField<T> * field, IBlock * const block ) >() ) :
      domain_decomposition::BlockDataCreator< FlagField<T> >( shared_ptr< field::DefaultBlockDataHandling< FlagField<T> > >(),
                                                              identifier, requiredSelectors, incompatibleSelectors )
   {
      if( alwaysInitialize )
      {
         auto dataHandling = make_shared< field::AlwaysInitializeBlockDataHandling< FlagField<T> > >( blocks, nrOfGhostLayers );
         dataHandling->addInitializationFunction( initFunction );
         this->dataHandling_ = dataHandling;
      }
      else
      {
         auto dataHandling = make_shared< field::DefaultBlockDataHandling< FlagField<T> > >( blocks, nrOfGhostLayers );
         dataHandling->addInitializationFunction( initFunction );
         this->dataHandling_ = dataHandling;
      }
   }
   
   Creator( const std::string & identifier = std::string(),
            const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet(),
            const uint_t nrOfGhostLayers = uint_t(1) ) :
      domain_decomposition::BlockDataCreator< FlagField<T> >( shared_ptr< DefaultBlockDataHandling< FlagField<T> > >(),
                                                              identifier, requiredSelectors, incompatibleSelectors )
   {
      static_assert( never_true<FlagField<T>>::value,
                     "This constructor of class 'field::Creator' is not supported anymore! Please use the other, new constructor." );
   }
};



} // namespace field
} // namespace walberla
