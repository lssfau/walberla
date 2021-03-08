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
//! \file FlagFieldMapping.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/FlagField.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"


namespace walberla {
namespace field {


template< typename FlagField_T, typename T >
class FlagFieldMapping : public vtk::BlockCellDataWriter<T,1>
{
private:
   using flag_t = typename FlagField_T::flag_t;
public:

   FlagFieldMapping( const ConstBlockDataID flagId, const std::string& id ) :
      vtk::BlockCellDataWriter<T,1>( id ), flagId_( flagId ), flagField_( NULL ) {}

   FlagFieldMapping( const ConstBlockDataID flagId, const std::string& id, const std::map< FlagUID, T > mapping ) :
      vtk::BlockCellDataWriter<T,1>( id ), flagId_( flagId ), flagField_( nullptr ), mapping_( mapping ) {}

   void addMapping( const FlagUID& flag, const T& value )
   {
      WALBERLA_ASSERT( mapping_.find( flag ) == mapping_.end() );
      mapping_[ flag ] = value;
   }

protected:

   void configure() override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( this->block_ );
      flagField_ = this->block_->template getData< FlagField_T >( flagId_ );

      for( auto mapping = mapping_.begin(); mapping != mapping_.end(); ++mapping )
      {
         if( flagField_->flagExists( mapping->first ) ) {
            flagMap_[ flagField_->getFlag( mapping->first ) ] = mapping->second;
         }
      }
   }

   T evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
      T result = 0;
      for( auto mapping = flagMap_.begin(); mapping != flagMap_.end(); ++mapping )
         if( flagField_->isFlagSet( x, y, z, mapping->first ) )
            result = static_cast<T>( result | mapping->second );
      return result;
   }

   const ConstBlockDataID flagId_;
   const FlagField_T*     flagField_;

   std::map< FlagUID, T > mapping_;
   std::map< flag_t , T > flagMap_;

}; // FlagFieldMapping


template<typename FieldType, typename TargetType=uint8_t>
class BinarizationFieldWriter : public vtk::BlockCellDataWriter<TargetType,1>
{
   using SrcType = typename FieldType::value_type;

public:
   BinarizationFieldWriter( const ConstBlockDataID fieldID, const std::string& id, SrcType mask) :
           vtk::BlockCellDataWriter<TargetType,1>( id ), fieldID_( fieldID ), field_( NULL ), mask_( mask ) {}

protected:

   void configure()  {
      WALBERLA_ASSERT_NOT_NULLPTR( this->block_ );
      field_ = this->block_->template getData< FieldType >( fieldID_ );
   }

   TargetType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
      if (field_->get(x,y,z) & mask_) {
         return TargetType(1);
      } else {
         return TargetType(0);
      }
   }
   const ConstBlockDataID fieldID_;
   const FieldType*       field_;

   SrcType mask_;
}; // BinaryFieldWriter



} // namespace field
} // namespace walberla
