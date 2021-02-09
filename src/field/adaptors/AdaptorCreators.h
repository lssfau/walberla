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
//! \file AdaptorCreators.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockDataHandling.h"

#include <functional>



namespace walberla {
namespace field {



template< typename Adaptor_T >
class AdaptorHandling : public blockforest::AlwaysInitializeBlockDataHandling< Adaptor_T >
{
public:

   typedef typename Adaptor_T::functor_t AdaptionFunction_T;

   AdaptorHandling( const ConstBlockDataID & adaptedFieldId, const AdaptionFunction_T & function ) :
      adaptedFieldId_( adaptedFieldId ), function_( function ) {}

   Adaptor_T * initialize( IBlock * const block ) override
   {
      typedef typename Adaptor_T::basefield_t AdaptedField_T;
      const AdaptedField_T * adaptedField = block->getData< AdaptedField_T >( adaptedFieldId_ );
      return new Adaptor_T( *adaptedField, function_ );
   }

private:

   ConstBlockDataID adaptedFieldId_;
   AdaptionFunction_T function_;
   
}; // class AdaptorHandling



template< typename Adaptor_T, typename BlockStorage_T > // BlockStorage_T will be deduced automatically
inline BlockDataID addFieldAdaptor( const shared_ptr< BlockStorage_T > & blocks,
                                    const ConstBlockDataID & adaptedFieldId,
                                    const std::string & identifier = std::string(),
                                    const typename Adaptor_T::functor_t & function = typename Adaptor_T::functor_t(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< AdaptorHandling<Adaptor_T> >( adaptedFieldId, function ), identifier, requiredSelectors, incompatibleSelectors );
}



template< typename Adaptor_T, typename BlockStorage_T > // BlockStorage_T will be deduced automatically
inline BlockDataID addFieldAdaptor( const shared_ptr< BlockStorage_T > & blocks,
                                    const ConstBlockDataID & adaptedFieldId,
                                    const std::string & identifier,
                                    const Set<SUID> & requiredSelectors,
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< AdaptorHandling<Adaptor_T> >( adaptedFieldId, typename Adaptor_T::functor_t() ), identifier, requiredSelectors, incompatibleSelectors );
}


/*// leads to problems when compiled with Intel compiler
template< typename Adaptor_T, typename BlockStorage_T > // BlockStorage_T will be deduced automatically
inline BlockDataID addFieldAdaptor( const BlockStorage_T & blocks,
                                    const ConstBlockDataID & adaptedFieldId,
                                    const std::string & identifier = std::string(),
                                    const typename Adaptor_T::functor_t & function = typename Adaptor_T::functor_t(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks.addBlockData( make_shared< AdaptorHandling<Adaptor_T> >( adaptedFieldId, function ), identifier, requiredSelectors, incompatibleSelectors );
}
*/


} // namespace field
} // namespace walberla
