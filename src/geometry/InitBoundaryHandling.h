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
//! \file Boundary.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief All boundary initializers combined
//
//======================================================================================================================

#pragma once


#include "initializer/BoundaryFromCellInterval.h"
#include "initializer/BoundaryFromDomainBorder.h"
#include "initializer/BoundaryFromImage.h"
#include "initializer/BoundaryFromVoxelFile.h"
#include "initializer/BoundaryFromBody.h"
#include "initializer/InitializationManager.h"

#include "structured/GrayScaleImage.h"
#include "structured/RGBAImage.h"

#include "core/config/Config.h"

namespace walberla {
namespace geometry {


//*******************************************************************************************************************
/*! Convenience function for setting up boundary handling via a block in the configuration file
*
*  This function groups together the most common initializers for a boundary handling.
*  For a detailed documentation on the format of the configuration block look at the initializer
*  documentation:
*
*  Initializer                                     | Functionality                                                                  | Subblock name
*  ------------------------------------------------|--------------------------------------------------------------------------------|----------------
*  initializer::BoundaryFromCellInterval           | \copybrief walberla::geometry::initializer::BoundaryFromCellInterval           | CellInterval
*  initializer::BoundaryFromDomainBorder           | \copybrief walberla::geometry::initializer::BoundaryFromDomainBorder           | Border
*  initializer::BoundaryFromImage<GrayScaleImage>  | \copybrief walberla::geometry::initializer::BoundaryFromImage<GrayScaleImage>  | GrayScaleImage
*  initializer::BoundaryFromImage<RGBAImage>       | \copybrief walberla::geometry::initializer::BoundaryFromImage<RGBAImage>       | RGBAImage
*  initializer::BoundaryFromVoxelFile              | \copybrief walberla::geometry::initializer::BoundaryFromVoxelFile              | VoxelFile
*  initializer::BoundaryFromBody                   | \copybrief walberla::geometry::initializer::BoundaryFromBody                   | Body
*
*  Example for a configuration block:
*
*  \code
*  Boundaries
*  {
*     Border { direction W;    walldistance -1;               Velocity0 {} }
*     Border { direction E;    walldistance -1;               Pressure1 {} }
*     Border { direction S,N;  walldistance -1;               NoSlip    {} }
*     Body   { shape Sphere;   radius 2; midpoint <10,10,10>; NoSlip    {} }
*
*     Image
*     {
*        file                 obstacle.png;
*        extrusionCoordinate  2;
*        rescaleToDomain      true;
*
*        BoundaryValueMapping {
*           NoSlip   {}
*           value    0;
*        }
*     }
*  }
*  \endcode
*/
//*******************************************************************************************************************
template<typename BoundaryHandling>
void initBoundaryHandling( StructuredBlockStorage & blocks, BlockDataID boundaryHandlingId,
                           const Config::BlockHandle & geometryBlock )
{
   using namespace geometry;
   using namespace initializer;

   using FromCellInterval = BoundaryFromCellInterval<BoundaryHandling>;
   using FromBorder = BoundaryFromDomainBorder<BoundaryHandling>;
   using FromVoxelFile = BoundaryFromVoxelFile<BoundaryHandling>;
   using FromBody = BoundaryFromBody<BoundaryHandling>;
   using FromGrayScaleImage = BoundaryFromImage<BoundaryHandling, GrayScaleImage>;
   using FromRGBAImage = BoundaryFromImage<BoundaryHandling, RGBAImage>;

   InitializationManager initManager( blocks.getBlockStorage() );


   initManager.registerInitializer( "CellInterval"  , shared_ptr<FromCellInterval   >(new FromCellInterval  ( blocks, boundaryHandlingId ) ));
   initManager.registerInitializer( "Border"        , shared_ptr<FromBorder         >(new FromBorder        ( blocks, boundaryHandlingId ) ));
   initManager.registerInitializer( "VoxelFile"     , shared_ptr<FromVoxelFile      >(new FromVoxelFile     ( blocks, boundaryHandlingId ) ));
   initManager.registerInitializer( "Body"          , shared_ptr<FromBody           >(new FromBody          ( blocks, boundaryHandlingId ) ));
   initManager.registerInitializer( "Image",          shared_ptr<FromGrayScaleImage >(new FromGrayScaleImage( blocks, boundaryHandlingId ) )); // for backwards compatibility
   initManager.registerInitializer( "GrayScaleImage", shared_ptr<FromGrayScaleImage >(new FromGrayScaleImage( blocks, boundaryHandlingId ) ));
   initManager.registerInitializer( "RGBAImage"     , shared_ptr<FromRGBAImage      >(new FromRGBAImage     ( blocks, boundaryHandlingId ) ));

   initManager.init( geometryBlock );
}

template<typename BoundaryHandling>
void initBoundaryHandling( StructuredBlockStorage & blocks, BlockDataID boundaryHandlingId, const Config & config )
{
   initBoundaryHandling<BoundaryHandling>( blocks, boundaryHandlingId, config.getBlock("Geometry") );
}



template<typename BoundaryHandling>
void setNonBoundaryCellsToDomain( StructuredBlockStorage & blocks, BlockDataID boundaryHandlingId )
{
   for( auto blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt )
   {
      auto handling = blockIt->template getData<BoundaryHandling>( boundaryHandlingId );
      auto flagField = handling->getFlagField();
      handling->fillWithDomain( flagField->nrOfGhostLayers() );
   }
}

template<typename FlagField_T>
void setNonBoundaryCellsToDomain( StructuredBlockStorage & blocks, BlockDataID flagFieldID,
                                  field::FlagUID fluidFlagID)
{
   for( auto blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt )
   {
      auto flagField = blockIt->template getData<FlagField_T>( flagFieldID );
      auto fluidFlag = flagField->getOrRegisterFlag(fluidFlagID);
      for( auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end(); ++it )
         if ( *it == 0 )
            addFlag(it, fluidFlag);
   }
}


} // namespace geometry
} // namespace walberla


