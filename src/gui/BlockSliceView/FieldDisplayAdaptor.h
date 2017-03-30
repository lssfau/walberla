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
//! \file FieldDisplayAdaptor.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DisplayAdaptor.h"
#include "field/GhostLayerField.h"


namespace walberla {
namespace gui {


   //*******************************************************************************************************************
   /*!
   * Abstract class, for displaying fields.
   *
   * Handles common functionality like fetching the field, generating slices etc.
   *
   *
   * field_t has to implement the concept of a GhostLayerField
   *
   * \ingroup gui
   */
   //*******************************************************************************************************************
   template<typename field_t>
   class FieldDisplayAdaptor : public DisplayAdaptor
   {
   public:
      FieldDisplayAdaptor( ConstBlockDataID fieldId )
         : field_(NULL), sliceDim_(-1), blockDataId_ ( fieldId )
      {}

      void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                      Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers );

   protected:
      Cell permuteCoordAccordingToSlice(  Cell input, int sliceCoord );
      void getFieldSizeInfo( CellInterval & ciOut, cell_idx_t & ghostLayerOut );

      const field_t *   field_;
      CellInterval      sliceInterval_;
      cell_idx_t        sliceDim_;
      ConstBlockDataID  blockDataId_;

   };

   template<typename field_t>
   void FieldDisplayAdaptor<field_t>::getFieldSizeInfo( CellInterval & ciOut, cell_idx_t & ghostLayerOut )
   {
      typedef GhostLayerField< typename field_t::value_type, field_t::F_SIZE > GhostField;

      const GhostField * glFieldPtr = dynamic_cast< const GhostField*> ( field_ );

      if ( glFieldPtr ) {
         ciOut = glFieldPtr->xyzSizeWithGhostLayer();
         ghostLayerOut = cell_idx_c( glFieldPtr->nrOfGhostLayers() );
      }
      else {
         ciOut = field_->xyzSize();
         ghostLayerOut = 0;
      }
   }




   template<typename field_t>
   void FieldDisplayAdaptor<field_t>::configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                                              Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers )
   {
      innerSize = Vector3<uint_t>(0,0,0);
      ghostLayers = 0;

      sliceDim_ = sliceDim;

      // get field pointer from BlockID, check if it is a ghost layer field
      field_ = block->getData<field_t>( blockDataId_ );

      // If field is not allocated on this block
      if ( ! field_ ) {
         sliceInterval_ = CellInterval();
         innerSize = Vector3<uint_t>(0,0,0);
         ghostLayers = 0;
         return;
      }

      CellInterval ci;
      getFieldSizeInfo( ci, ghostLayers  );

      // if slice value is not in the valid range nothing is displayed
      if ( sliceValue < ci.min()[ uint_c(sliceDim)] || sliceValue > ci.max()[uint_c(sliceDim)] ) {
         sliceInterval_ = CellInterval();
         return;
      }

      // Determine size
      Cell cSize = Cell ( field_->xSize(), field_->ySize(), field_->zSize() );
      cSize = permuteCoordAccordingToSlice( cSize, sliceDim );
      innerSize = Vector3<uint_t> ( uint_c( cSize[0u] ), uint_c( cSize[1u] ), uint_c( cSize[2u] ) );

      // create slice interval
      sliceInterval_ = ci;
      sliceInterval_.min()[ uint_c(sliceDim) ] = sliceValue;
      sliceInterval_.max()[ uint_c(sliceDim) ] = sliceValue;
   }


   template<typename field_t>
   Cell FieldDisplayAdaptor<field_t>::permuteCoordAccordingToSlice(  Cell input, int sliceCoord )
   {
      static const int permutation[3][3] = { {1,2,0}, {0,2,1}, {0,1,2} };

      Cell coord;
      for( uint_t i=0; i<3; ++i )
         coord[i] = input[ uint_c( permutation[sliceCoord][i] ) ];

      return coord;
   }




} // namespace gui
} // namespace walberla


