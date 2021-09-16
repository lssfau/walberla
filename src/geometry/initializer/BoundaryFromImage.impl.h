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
//! \file BoundaryFromImage.impl.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#include "core/Abort.h"
#include "core/math/Vector3.h"
#include "field/GhostLayerField.h"

#include <limits>
#include "BoundaryFromImage.h"


namespace walberla {
namespace geometry {
namespace initializer {


   //===================================================================================================================
   //
   //  Helper functions
   //
   //===================================================================================================================

   static inline void transform3Dto2D( Cell & cell, uint_t extrusionCoord )  {
      switch ( extrusionCoord ) {
         case 0:   cell = Cell( cell[1], cell[2], cell[0] );   break;
         case 1:   cell = Cell( cell[0], cell[2], cell[1] );   break;
         case 2:                                               break;
         default:  WALBERLA_ASSERT( false);                             break;
      }
   }

   static inline void transform2Dto3D( Cell & cell, uint_t extrusionCoord )  {
      switch ( extrusionCoord ) {
         case 0:   cell = Cell( cell[2], cell[0], cell[1] );   break;
         case 1:   cell = Cell( cell[0], cell[2], cell[1] );   break;
         case 2:                                               break;
         default:  WALBERLA_ASSERT( false);                             break;
      }
   }

   static inline void transform2Dto3D( CellInterval & ci, uint_t extrusionCoord )
   {
      transform2Dto3D( ci.min(), extrusionCoord );
      transform2Dto3D( ci.max(), extrusionCoord );
   }


   //===================================================================================================================
   //
   //  BoundaryFromImage
   //
   //===================================================================================================================

   template<typename Handling, typename Image_T >
   BoundaryFromImage<Handling, Image_T>::BoundaryFromImage( StructuredBlockStorage & structuredBlockStorage,
                                                            BlockDataID handlerBlockDataID )
      : structuredBlockStorage_( structuredBlockStorage ),
        handlingID_( handlerBlockDataID )
   {}


   template<typename Handling, typename Image_T >
   void BoundaryFromImage<Handling, Image_T>::init( const Image_T & img, uint_t extrusionCoord,
                                                    BoundarySetters & boundarySetters,
                                                    cell_idx_t xOffset, cell_idx_t yOffset,
                                                    cell_idx_t lowerExtrusionLimit,
                                                    cell_idx_t upperExtrusionLimit )
   {
      if ( boundarySetters.empty() )
         return;

      // the extrusion limit can not be left at numeric_limits::min/max since
      // the blockStorage.transformGlobalToLocal would overflow
      cell_idx_t maxNumberOfCells = cell_idx_c(
                                    std::max( structuredBlockStorage_.getNumberOfXCells(),
                                    std::max( structuredBlockStorage_.getNumberOfYCells(),
                                              structuredBlockStorage_.getNumberOfZCells() ) ) );

      if ( upperExtrusionLimit > maxNumberOfCells )
         upperExtrusionLimit = maxNumberOfCells;


      for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
      {
         for( auto setterIt = boundarySetters.begin(); setterIt != boundarySetters.end(); ++setterIt )
            setterIt->second.configure( *blockIt, handlingID_  );

         const auto flagField = BoundarySetter<Handling>::getFlagField(*blockIt, handlingID_);

         const cell_idx_t width  = cell_idx_c( img.width() ) ;
         const cell_idx_t height = cell_idx_c( img.height()) ;
         CellInterval imgInterval( xOffset,          yOffset,            lowerExtrusionLimit,
                                   xOffset+ width-1, yOffset + height-1, upperExtrusionLimit-1 );

         const cell_idx_t ghostLayers = cell_idx_c( flagField->nrOfGhostLayers() );

         if ( lowerExtrusionLimit < -ghostLayers )
            lowerExtrusionLimit = -ghostLayers;

         transform2Dto3D( imgInterval, extrusionCoord );

         structuredBlockStorage_.transformGlobalToBlockLocalCellInterval( imgInterval, *blockIt );

         imgInterval.intersect( flagField->xyzSizeWithGhostLayer() );
         for( auto cellIt = flagField->beginSliceXYZ( imgInterval ); cellIt != flagField->end(); ++cellIt )
         {
            Cell posInImage = cellIt.cell();
            structuredBlockStorage_.transformBlockLocalToGlobalCell( posInImage, *blockIt );
            transform3Dto2D( posInImage, extrusionCoord );

            const Cell posInField = cellIt.cell();

            const pixel_t imgValue = img.getElement( posInImage[0] - xOffset , posInImage[1] - yOffset );

            auto boundarySetter = boundarySetters.find( imgValue );
            if ( boundarySetter != boundarySetters.end() ) {
               boundarySetter->second.set( posInField[0], posInField[1], posInField[2] );
            }
         }
      }
   }


   template<typename Handling, typename Image_T >
   void BoundaryFromImage<Handling, Image_T>::init ( const Config::BlockHandle & block )
   {
      // Read Boundary value mapping
      if ( structuredBlockStorage_.begin() == structuredBlockStorage_.end() )
         return; // No blocks on this process

      if ( ! block.isDefined( "file") )
         WALBERLA_ABORT( "Missing Parameter 'file' in scalar field to image block" );

      const std::string & file = block.getParameter< std::string > ( "file" );

      uint_t extrusionCoordinate = block.getParameter<uint_t>( "extrusionCoordinate", 2 );

      const cell_idx_t minCellIdx = std::numeric_limits<cell_idx_t>::min();
      const cell_idx_t maxCellIdx = std::numeric_limits<cell_idx_t>::max();

      cell_idx_t lowerExtrusionLimit = block.getParameter<cell_idx_t>( "lowerExtrusionLimit", minCellIdx );
      cell_idx_t upperExtrusionLimit = block.getParameter<cell_idx_t>( "upperExtrusionLimit", maxCellIdx );

      cell_idx_t xOffset = block.getParameter<cell_idx_t>( "xOffset", 0 );
      cell_idx_t yOffset = block.getParameter<cell_idx_t>( "yOffset", 0 );


      Config::Blocks mappingBlocks;
      block.getBlocks( "BoundaryValueMapping", mappingBlocks, 1 );

      BoundarySetters mapping;

      for( auto mappingBlock = mappingBlocks.begin(); mappingBlock != mappingBlocks.end(); ++mappingBlock )
      {
         const std::string pixelValueStr =  mappingBlock->template getParameter<std::string> ( "value" );
         auto pixelValue = Image_T::pixelValueFromString( pixelValueStr );

         if( mapping.find( pixelValue )  != mapping.end() )
            WALBERLA_ABORT( "BoundaryFromImage: Duplicate BoundaryValueMapping for pixel value " << pixelValue );

         BoundarySetter<Handling> boundarySetter;
         boundarySetter.setConfigBlock( *mappingBlock );

         mapping[pixelValue] = boundarySetter;
      }

      Image_T img ( file );
      if ( block.isDefined( "rescaleToDomain") )
         img = resizeImageToDomain( img, extrusionCoordinate );

      init( img, extrusionCoordinate, mapping,  xOffset, yOffset, lowerExtrusionLimit, upperExtrusionLimit );
   }

   template<typename Handling, typename Image_T >
   Image_T BoundaryFromImage<Handling, Image_T>::resizeImageToDomain( const Image_T & img, uint_t extrusionCoord ) const
   {
      Cell domainSize ( structuredBlockStorage_.getNumberOfXCells(),
                        structuredBlockStorage_.getNumberOfYCells(),
                        structuredBlockStorage_.getNumberOfZCells() );

      transform3Dto2D( domainSize, extrusionCoord );
      return img.getResizedImage( uint_c(domainSize[0]), uint_c(domainSize[1]), false );
   }



} // namespace initializer
} // namespace geometry
} // namespace walberla


