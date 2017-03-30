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
//! \file ScalarFieldFromGrayScaleImage.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "ScalarFieldFromGrayScaleImage.h"
#include "core/Abort.h"
#include "core/math/Vector3.h"
#include "field/GhostLayerField.h"

#include <limits>


namespace walberla {
namespace geometry {
namespace initializer {


   //===================================================================================================================
   //
   //  Helper functions
   //
   //===================================================================================================================

   static void transform3Dto2D( Cell & cell, uint_t extrusionCoord )  {
      switch ( extrusionCoord ) {
         case 0:   cell = Cell( cell[1], cell[2], cell[0] );   break;
         case 1:   cell = Cell( cell[0], cell[2], cell[1] );   break;
         case 2:                                               break;
         default:  WALBERLA_ASSERT( false);                             break;
      }
   }

   static void transform2Dto3D( Cell & cell, uint_t extrusionCoord )  {
      switch ( extrusionCoord ) {
         case 0:   cell = Cell( cell[2], cell[0], cell[1] );   break;
         case 1:   cell = Cell( cell[0], cell[2], cell[1] );   break;
         case 2:                                               break;
         default:  WALBERLA_ASSERT( false);                             break;
      }
   }

   static void transform2Dto3D( CellInterval & ci, uint_t extrusionCoord )
   {
      transform2Dto3D( ci.min(), extrusionCoord );
      transform2Dto3D( ci.max(), extrusionCoord );
   }


   //===================================================================================================================
   //
   //  ScalarFieldFromGrayScaleImage
   //
   //===================================================================================================================


   ScalarFieldFromGrayScaleImage::ScalarFieldFromGrayScaleImage( StructuredBlockStorage & structuredBlockStorage,
                                                                 BlockDataID scalarFieldID )
      : structuredBlockStorage_( structuredBlockStorage ),
        scalarFieldID_( scalarFieldID )
   {}


   void ScalarFieldFromGrayScaleImage::init ( BlockStorage & , const Config::BlockHandle & block )
   {
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


      GrayScaleImage img ( file );
      if ( block.isDefined( "scaleToDomain") )
      {
         if ( xOffset != 0 || yOffset != 0 )
            WALBERLA_ABORT("Error when loading image " << file
                           << " : When 'scaleToDomain' is specified [xy]Offset has to be zero" );

         init( img, extrusionCoordinate, true, lowerExtrusionLimit, upperExtrusionLimit );
      }
      else
      {
         init( img, extrusionCoordinate, xOffset, yOffset, lowerExtrusionLimit, upperExtrusionLimit );
      }
   }


   void ScalarFieldFromGrayScaleImage::init( const GrayScaleImage & img, uint_t extrusionCoord,
                                             cell_idx_t xOffset, cell_idx_t yOffset,
                                             cell_idx_t lowerExtrusionLimit,
                                             cell_idx_t upperExtrusionLimit )
   {

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
         GhostLayerField<real_t,1> * f = blockIt->getData<GhostLayerField<real_t,1> >( scalarFieldID_ );

         const cell_idx_t width  = cell_idx_c( img.width() ) ;
         const cell_idx_t height = cell_idx_c( img.height()) ;
         CellInterval imgInterval( xOffset,          yOffset,            lowerExtrusionLimit,
                                   xOffset+ width-1, yOffset + height-1, upperExtrusionLimit-1 );


         cell_idx_t lowestGlobalCoord = - cell_idx_c( f->nrOfGhostLayers() );
         if ( lowerExtrusionLimit < lowestGlobalCoord )
            lowerExtrusionLimit = lowestGlobalCoord;


         transform2Dto3D( imgInterval, extrusionCoord );

         structuredBlockStorage_.transformGlobalToBlockLocalCellInterval( imgInterval, *blockIt );

         imgInterval.intersect( f->xyzSizeWithGhostLayer() );

         for( auto cellIt = f->beginSliceXYZ( imgInterval ); cellIt != f->end(); ++cellIt )
         {
            Cell curPos = cellIt.cell();
            structuredBlockStorage_.transformBlockLocalToGlobalCell( curPos, *blockIt );
            transform3Dto2D( curPos, extrusionCoord );
            *cellIt = img( curPos[0] - xOffset , curPos[1] - yOffset );
         }
      }
   }

   void ScalarFieldFromGrayScaleImage::init( const GrayScaleImage & img, uint_t extrusionCoord, bool rescaleToDomain,
                                             cell_idx_t lowerExtrusionLimit, cell_idx_t upperExtrusionLimit )
   {
      if ( rescaleToDomain )
      {
         Cell domainSize ( structuredBlockStorage_.getNumberOfXCells(),
                           structuredBlockStorage_.getNumberOfYCells(),
                           structuredBlockStorage_.getNumberOfZCells() );

         transform3Dto2D( domainSize, extrusionCoord );

         init(  img.getResizedImage( uint_c( domainSize[0] ), uint_c( domainSize[1] ) ), extrusionCoord, 0,0, lowerExtrusionLimit, upperExtrusionLimit );
      }
      else
         init( img, extrusionCoord, 0, 0, lowerExtrusionLimit, upperExtrusionLimit );
   }



} // namespace initializer
} // namespace geometry
} // namespace walberla


