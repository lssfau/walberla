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
//! \file FlagFieldCellFilter.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/FlagField.h"
#include "core/cell/CellSet.h"
#include "core/debug/Debug.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace field {



template< typename FlagField_T >
class FlagFieldCellFilter {
private:
   using flag_t = typename FlagField_T::flag_t;
public:

   FlagFieldCellFilter( const ConstBlockDataID flags ) : flagField_( flags ) {}

   FlagFieldCellFilter( const ConstBlockDataID flags, const FlagUID & filteredFlag )
      : flagField_( flags ), filteredFlags_( 1, filteredFlag ) {}

   void addFlag( const FlagUID& flag ) { filteredFlags_.push_back( flag ); }

   void operator()( CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage, const uint_t ghostLayers = uint_t(0) ) const
   {
      const FlagField_T* flagField = block.getData< FlagField_T >( flagField_ );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      std::vector< FlagUID > existingFlags;
      for( auto flag = filteredFlags_.begin(); flag != filteredFlags_.end(); ++flag )
         if( flagField->flagExists( *flag ) )
            existingFlags.push_back( *flag );

      if( !existingFlags.empty() )
      {
         flag_t filterMask( flag_t(0) );
         for( auto flag = existingFlags.begin(); flag != existingFlags.end(); ++flag )
            filterMask = static_cast< flag_t >( filterMask | flagField->getFlag( *flag ) );

         const cell_idx_t gl    = cell_idx_c( ghostLayers );
         const cell_idx_t begin = cell_idx_c( -1 ) * gl;

         for( cell_idx_t z = begin; z < cell_idx_c( storage.getNumberOfZCells(block) ) + gl; ++z )
            for( cell_idx_t y = begin; y < cell_idx_c( storage.getNumberOfYCells(block) ) + gl; ++y )
               for( cell_idx_t x = begin; x < cell_idx_c( storage.getNumberOfXCells(block) ) + gl; ++x )
                  if( flagField->isPartOfMaskSet( x, y, z, filterMask ) ) filteredCells.insert( x, y, z );
      }
   }

private:

   const ConstBlockDataID flagField_;
   std::vector< FlagUID > filteredFlags_;

}; // class FlagFieldCellFilter



} // namespace field
} // namespace walberla
