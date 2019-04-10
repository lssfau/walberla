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
//! \file OverlapFieldFromBody.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "OverlapFieldFromBody.h"

#include "geometry/bodies/AABBBody.h"
#include "geometry/bodies/BodyFromConfig.h"
#include "geometry/bodies/BodyOverlapFunctions.h"
#include "geometry/bodies/Cylinder.h"
#include "geometry/bodies/Ellipsoid.h"
#include "geometry/bodies/Sphere.h"
#include "geometry/bodies/Torus.h"

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/StringUtility.h"
#include "field/GhostLayerField.h"


namespace walberla {
namespace geometry {
namespace initializer {


   OverlapFieldFromBody::OverlapFieldFromBody( StructuredBlockStorage & structuredBlockStorage, BlockDataID scalarFieldID,
                                             const std::string & addKeyword, const std::string & subtractKeyword )
      : structuredBlockStorage_( structuredBlockStorage ),
        scalarFieldID_( scalarFieldID ),
        addKeyword_ ( addKeyword ),
        subtractKeyword_ ( subtractKeyword )
   {
   }


   void OverlapFieldFromBody::init( BlockStorage & /*blockStorage*/, const Config::BlockHandle & blockHandle )
   {
      Config::Blocks subBlocks;
      blockHandle.getBlocks( subBlocks );

      if ( blockHandle.isDefined("initialFill") ) {
         std::string initialFill = blockHandle.getParameter<std::string>("initialFill");
         if ( initialFill == addKeyword_ )
         {
            for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
               blockIt->getData<GhostLayerField<real_t,1> > ( scalarFieldID_ )->setWithGhostLayer( real_t(1.0) );
         }
         else if ( initialFill == subtractKeyword_ ) {
            for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
               blockIt->getData<GhostLayerField<real_t,1> > ( scalarFieldID_ )->setWithGhostLayer( real_t(0.0) );
         }
         else {
            WALBERLA_ABORT("Unknown value of initialFill. Valid values are " << addKeyword_ << "," << subtractKeyword_ );
         }
      }

      for ( auto it = subBlocks.begin(); it != subBlocks.end(); ++it )
      {
         Config::BlockHandle subBlock = *it;
         bool addOrSubtract = true;

         bool addDefined      = subBlock.isDefined( addKeyword_ );
         bool subtractDefined = subBlock.isDefined( subtractKeyword_ );

         if ( addDefined && subtractDefined )
            WALBERLA_ABORT( "Specify only one of " << addKeyword_ << " and " << subtractKeyword_ << "!\n"
                            << "Both are defined in " << blockHandle.getKey() );

         if ( subtractDefined )
            addOrSubtract = false;

         std::string shape = subBlock.getParameter<std::string>("shape");

         uint_t superSamplingDepth = subBlock.getParameter<uint_t>("superSamplingDepth", uint_t(4) );

         if      ( string_icompare( shape, "Sphere"   ) == 0 )  init ( sphereFromConfig   ( subBlock ), addOrSubtract, superSamplingDepth );
         else if ( string_icompare( shape, "Cylinder" ) == 0 )  init ( cylinderFromConfig ( subBlock ), addOrSubtract, superSamplingDepth );
         else if ( string_icompare( shape, "Torus"    ) == 0 )  init ( torusFromConfig    ( subBlock ), addOrSubtract, superSamplingDepth );
         else if ( string_icompare( shape, "Ellipsoid") == 0 )  init ( ellipsoidFromConfig( subBlock ), addOrSubtract, superSamplingDepth );
         else if ( string_icompare( shape, "Box"      ) == 0 )  init ( AABBFromConfig     ( subBlock ), addOrSubtract, superSamplingDepth );
         else
         {
            WALBERLA_ABORT( "Unknown Block " << subBlock.getKey() << " in block " << blockHandle.getKey() << "\n"
                            << "Allowed blocks are 'Sphere', 'Ellipsoid' and 'Box'" );
         }

      }
   }

} // namespace initializer
} // namespace geometry
} // namespace walberla
