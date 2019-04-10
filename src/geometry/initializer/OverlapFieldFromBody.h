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
//! \file OverlapFieldFromBody.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Initializer.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "geometry/bodies/BodyFromConfig.h"
#include "geometry/bodies/BodyOverlapFunctions.h"

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "field/GhostLayerField.h"


namespace walberla {
namespace geometry {
namespace initializer {



   //*******************************************************************************************************************
   /*! Initializes a scalar field from a geometric body
   *
   * Currently supported are Sphere, Ellipsoid and Box (= AABB)
   *
   * Example:
   * \verbatim
          <InitializerUID> {
               initialFill : drop;
               someArbitraryId {
                  shape sphere;
                  bubble;

                  midpoint < 3,4,5>;
                  radius 4;
               }
               object2 {
                  shape box;
                  drop;

                  min <1,2,3>;
                  max <3,4,5>;
               }
               object3_ellipse {
                  bubble;
                  shape ellipsoid;
                  midpoint < 3,4,2>;
                  axis1    <1,0,0>;
                  axis2    <0,1,0>;
                  radii    <1,1,4>;
               }
          }
   * \endverbatim
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   class OverlapFieldFromBody : public Initializer
   {
   public:

      /*************************************************************************************************************//**
      * Constructor
      *
      * \param scalarFieldID    the scalar field to initialize
      * \param addKeyword       used when parsing a configuration block, and determining the addOrSubtract parameter for
      *                         ScalarFieldFromBody::init(const Body&, bool) -> see documentation of this function.
      *                         If the addKeyword is defined in the block, the overlapFraction is added to the cells
      *                         otherwise subtracted
      * \param subtractKeyword  see above
      *
      *****************************************************************************************************************/
      OverlapFieldFromBody( StructuredBlockStorage & structuredBlockStorage, BlockDataID scalarFieldID,
                           const std::string & addKeyword = "drop", const std::string & subtractKeyword = "bubble" );



      /*************************************************************************************************************//**
      * Interface implementation for Initializer - sets a body on a scalar field with options from configuration file
      *
      *****************************************************************************************************************/
      virtual void init( BlockStorage & blockStorage, const Config::BlockHandle & blockHandle );



      /*************************************************************************************************************//**
      * Sets a body on the scalar field
      *
      * \param body The body object - has to implement either overlapFraction(...), or contains(...)
      *             see BodyOverlapFunctions for detailed body concept
      * \param addOrSubtract if true the overlap between body and cell (fraction between 0 and 1) is
      *                      added to scalar field, and if result is greater 1, then the value is set to one
      *                      if false, the overlap is subtracted, and minimum value is zero

      *  Supported bodies are Sphere, Ellipsoid, AABB.
      *  To add a new supported body implement concept defined in BodyOverlapFunctions.h, and
      *  add an explicit template instantiation in ScalarFieldFromBody.cpp for the new body.
      *
      *****************************************************************************************************************/
      template<typename Body>
      void init( const Body & body, bool addOrSubtract, uint_t superSamplingDepth=4 );


   protected:

      StructuredBlockStorage & structuredBlockStorage_;
      BlockDataID              scalarFieldID_;

      std::string              addKeyword_;
      std::string              subtractKeyword_;
   };

   template< typename Body >
   void OverlapFieldFromBody::init( const Body & body, bool addOrSubtract, uint_t superSamplingDepth )
   {
      const real_t dx = structuredBlockStorage_.dx();
      const real_t dy = structuredBlockStorage_.dy();
      const real_t dz = structuredBlockStorage_.dz();
      const Vector3<real_t> dxVec(dx, dy, dz);

      for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
      {
         IBlock * block = &(*blockIt);

         GhostLayerField<real_t,1> * ff = block->getData<GhostLayerField<real_t,1> >( scalarFieldID_ );
         auto gl = cell_idx_c( ff->nrOfGhostLayers() );

         // If Block (extended with ghost layers) does not intersect body - skip the complete block
         AABB blockBB = block->getAABB();
         blockBB.extend( math::Vector3< real_t >( dx * real_c( gl ), dy * real_c( gl ), dz * real_c( gl ) ) );
         if( fastOverlapCheck( body, blockBB ) == geometry::COMPLETELY_OUTSIDE )
            continue;

         AABB firstCellBB;
         structuredBlockStorage_.getBlockLocalCellAABB( *block, ff->beginWithGhostLayer().cell(), firstCellBB );
         Vector3<real_t> firstCellMidpoint;
         for( uint_t i = 0; i < 3; ++i )
            firstCellMidpoint[i] = firstCellBB.min(i) + real_t(0.5) * firstCellBB.size(i);

         Vector3<real_t> currentMidpoint;
         currentMidpoint[2] = firstCellMidpoint[2];
         for ( cell_idx_t z = -gl; z < cell_idx_c(ff->zSize())+gl; ++z, currentMidpoint[2] += dz )
         {
            currentMidpoint[1] = firstCellMidpoint[1];
            for ( cell_idx_t y = -gl; y < cell_idx_c(ff->ySize())+gl; ++y, currentMidpoint[1] += dy )
            {
               currentMidpoint[0] = firstCellMidpoint[0];
               for( cell_idx_t x = -gl; x < cell_idx_c(ff->xSize())+gl; ++x, currentMidpoint[0] += dx )
               {
                  real_t overlap = overlapFraction( body, currentMidpoint, dxVec, superSamplingDepth );
                  real_t & val = ff->get(x,y,z);
                  WALBERLA_ASSERT( val >=0 && val <= 1);

                  if ( addOrSubtract ) {
                     val += overlap;
                     val = std::min( val, real_t(1) );
                  }
                  else {
                     val -= overlap;
                     val = std::max( val, real_t(0) );
                  }

               }
            }
         }
      }
   }




} // namespace initializer
} // namespace geometry
} // namespace walberla



