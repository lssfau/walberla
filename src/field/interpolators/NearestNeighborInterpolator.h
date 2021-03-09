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
//! \file NearestNeighborInterpolator.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace field {


   //*******************************************************************************************************************
   /*!
   * Nearest Neighbor interpolation for GhostLayerFields
   *
   * \ingroup field
   *
   * Provides access to a field with real valued coordinates. The value of the nearest neighbor is
   * returned.
   *
   * Template parameter GlField has to be a GhostLayerField or a class that implements the methods
   *    CellInterval xyzSizeWithGhostLayer()
   *    ResType      get( cell_idx_t x, cell_idx_t y, cell_idx_t z )
   *
   * Coordinate System:
   *     \image html field/doc/interpolationCoordinates.png
   *
   */
   //*******************************************************************************************************************
   template< typename GlField >
   class NearestNeighborInterpolator
   {
   public:

      static const uint_t F_SIZE = GlField::F_SIZE;
      using value_type = typename GlField::value_type;

      NearestNeighborInterpolator( const GlField & fieldToInterpolate )
         : field_( fieldToInterpolate )
      {
         const CellInterval size = field_.xyzSizeWithGhostLayer();

         validRegion_ = AABB( real_c(size.xMin()    ), real_c(size.yMin()    ), real_c(size.zMin()    ),
                              real_c(size.xMax() + 1), real_c(size.yMax() + 1), real_c(size.zMax() + 1) );

      }

      inline value_type operator()( real_t x, real_t y, real_t z )
      {
         WALBERLA_ASSERT( F_SIZE == uint_t(1u) );
         return (*this)(x,y,z,0);
      }

      inline value_type operator()( real_t x, real_t y, real_t z, cell_idx_t f )
      {
         WALBERLA_ASSERT( validRegion_.contains( x,y,z ) );

         const cell_idx_t ciX = x < real_t(0) ? cell_idx_c( x - real_t(1) ) : cell_idx_c( x );
         const cell_idx_t ciY = y < real_t(0) ? cell_idx_c( y - real_t(1) ) : cell_idx_c( y );
         const cell_idx_t ciZ = z < real_t(0) ? cell_idx_c( z - real_t(1) ) : cell_idx_c( z );

         return field_.get( ciX, ciY, ciZ, f );
      }

      const AABB & getValidRegion() const {
         return validRegion_;
      }

   protected:
      const GlField & field_;
      AABB validRegion_;
   };


} // namespace field
} // namespace walberla




