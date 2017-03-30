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
//! \file FieldClone.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Class holds one copy for each size of the given field
//
//======================================================================================================================

#pragma once

#include "SwapableCompare.h"
#include "core/debug/Debug.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include <map>
#include <set>


namespace walberla {
namespace field {


//*******************************************************************************************************************
/*! Helper class for the common scenario where a 'scratch field' is needed
*
* Consider for example a LBM stream operation: Usually two fields are used, one source field
* to read the values and one destination field where the result is written to.
* In the end both fields are swapped.
* Now if one has multiple fields per block one scratch field is enough if all fields have
* the same size.
*
* This class helps in this scenario by providing that extra scratch field.
* The first time get() is called, a clone of the field is created and returned.
* Further calls return the same scratch field.
*
* If fields have different sizes one scratch field is held for each size.
*
* \tparam WorkField  field that has to be derived from GhostLayerField
* \tparam constSize  have all fields the same size?
*/
//*******************************************************************************************************************
template< typename WorkField, bool constSize = false >
class FieldClone
{

public:
   FieldClone( ConstBlockDataID fieldID ) : fieldID_( fieldID ), dstField_( NULL ){}

   ~FieldClone(){
      // Free allocated temporary fields
      for ( auto dst = dstFieldSet_.begin(); dst != dstFieldSet_.end(); ++dst )
         delete *dst;

      if( constSize && dstField_ )
         delete dstField_;
   }

   inline WorkField* get( IBlock* block )
   {
      bool created;
      return get( block, created );
   }

   inline WorkField* get( IBlock* block, bool& created )
   {
      created = false;

      WorkField* src = block->getData<WorkField> ( fieldID_ );

      if( constSize )
      {
         if( !dstField_ ){
            dstField_ = src->cloneUninitialized();
            created = true;
         }
         WALBERLA_ASSERT( src->hasSameAllocSize( *dstField_ ) );
         return dstField_;
      }
      else
      {
         // Search for an already allocated destination field of the same size as our given src field
         auto dst = dstFieldSet_.find ( src );
         if( dst != dstFieldSet_.end() ) {
            WALBERLA_ASSERT( src->hasSameAllocSize( **dst ) );
            return *dst;
         }

         // If none was found, allocate a new one
         WorkField* newDst = src->cloneUninitialized( );
         WALBERLA_ASSERT_NOT_NULLPTR( newDst );
         dstFieldSet_.insert( newDst );
         created = true;

         return newDst;
      }
    }

private:
   const BlockDataID fieldID_;   ///< field which is processed

   WorkField*                                           dstField_;      ///< field pointer for constant size
   std::set< WorkField*, SwapableCompare<WorkField*> >  dstFieldSet_;   ///< dst field for every field-size occurring on this process
};



/*
template< typename WorkField, bool constSize = false >
class FieldCreator
{

public:
   FieldCreator() : dstField_( NULL ){}

   ~FieldCreator(){
      // Free allocated temporary fields
      for ( auto dst = dstFieldMap_.begin(); dst != dstFieldMap_.end(); ++dst )
         delete dst->second;

      if( constSize && dstField_ )
         delete dstField_;
   }

   inline WorkField* get( IBlock* block, uint_t xSize, uint_t ySize, uint_t zSize )
   {
      bool created;
      return get( block, xSize, ySize, zSize, created );
   }

   inline WorkField* get( IBlock* block, uint_t xSize, uint_t ySize, uint_t zSize, bool& created )
   {
      created = false;

      if( constSize )
      {
         if( !dstField_ ){
            WALBERLA_ASSERT( xSize>0u && ySize>0u && zSize>0u );
            dstField_ = new WorkField( xSize, ySize, zSize );
            created = true;
         }
         return dstField_;
      }
      else
      {
         // Search for an already allocated destination field of the same size as our given src field
         auto dst = dstFieldMap_.find ( block );
         if( dst != dstFieldMap_.end() ) {
            return dst->second;
         }

         // If none was found, allocate a new one
         WALBERLA_ASSERT( xSize>0u && ySize>0u && zSize>0u );
         WorkField* newDst = new WorkField( xSize, ySize, zSize );
         WALBERLA_ASSERT_NOT_NULLPTR( newDst );
         dstFieldMap_[block] = newDst;
         created = true;

         return newDst;
      }
   }

private:
   WorkField*                       dstField_;      ///< field pointer for constant size
   std::map< IBlock*, WorkField* >  dstFieldMap_;   ///< dst field for every field-size occurring on this process
};
*/

} // namespace field
} // namespace walberla

