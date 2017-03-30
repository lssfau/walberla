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
//! \file TrilinearInterpolator.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "TrilinearInterpolatorFwd.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/GhostLayerField.h"


namespace walberla {
namespace field {

   //*******************************************************************************************************************
   /*!
   * Trilinear interpolation for GhostLayerFields
   *
   * \ingroup field
   *
   * Template parameter GlField has to be a GhostLayerField or a class that implements the methods
   *     CellInterval         xyzSizeWithGhostLayer()
   *     GlField::value_type  get( cell_idx_t x, cell_idx_t y, cell_idx_t z )
   *
   *
   * Provides access to a field with real valued coordinates.
   * Values are interpolated with a trilinear scheme.
   *
   * Coordinate System:
   *     \image html field/doc/interpolationCoordinates.png
   */
   //*******************************************************************************************************************
   template< typename GlField, typename ResType /*= real_t*/, typename FactorType /*= real_t*/ >
   class TrilinearInterpolator
   {
   public:

      static const uint_t F_SIZE = GlField::F_SIZE;


      TrilinearInterpolator( const GlField & fieldToInterpolate )
         : field_( fieldToInterpolate )
      {
         const CellInterval size = field_.xyzSizeWithGhostLayer();

         validRegion_ = AABB( real_c(size.xMin())             + real_c(0.5), real_c(size.yMin())             + real_c(0.5), real_c(size.zMin())             + real_c(0.5),
                              real_c(size.xMax()) + real_t(1) - real_c(0.5), real_c(size.yMax()) + real_t(1) - real_c(0.5), real_c(size.zMax()) + real_t(1) - real_c(0.5) );
      }

      inline ResType operator()( real_t x, real_t y, real_t z ) const
      {
         WALBERLA_ASSERT( F_SIZE == uint_t(1u) );
         return (*this)(x,y,z,0);
      }


      inline ResType operator()( real_t x, real_t y, real_t z, cell_idx_t f ) const
      {

         WALBERLA_ASSERT( validRegion_.contains( x,y,z ), "Point (" << x << "," << y << "," << z << ") not in valid region " << getValidRegion()  );

         const cell_idx_t xMin = x < real_c(0.5) ? cell_idx_c( x - real_c(1.5) ) : cell_idx_c ( x - real_c(0.5) );
         const cell_idx_t yMin = y < real_c(0.5) ? cell_idx_c( y - real_c(1.5) ) : cell_idx_c ( y - real_c(0.5) );
         const cell_idx_t zMin = z < real_c(0.5) ? cell_idx_c( z - real_c(1.5) ) : cell_idx_c ( z - real_c(0.5) );

         const FactorType xDiff = numeric_cast< FactorType >( x - ( real_c(xMin) + real_c(0.5) ) );
         const FactorType yDiff = numeric_cast< FactorType >( y - ( real_c(yMin) + real_c(0.5) ) );
         const FactorType zDiff = numeric_cast< FactorType >( z - ( real_c(zMin) + real_c(0.5) ) );

         WALBERLA_ASSERT_GREATER_EQUAL( xDiff, numeric_cast<FactorType>(0) );
         WALBERLA_ASSERT_GREATER_EQUAL( yDiff, numeric_cast<FactorType>(0) );
         WALBERLA_ASSERT_GREATER_EQUAL( zDiff, numeric_cast<FactorType>(0) );

         WALBERLA_ASSERT_LESS_EQUAL( xDiff, numeric_cast<FactorType>(1) );
         WALBERLA_ASSERT_LESS_EQUAL( yDiff, numeric_cast<FactorType>(1) );
         WALBERLA_ASSERT_LESS_EQUAL( zDiff, numeric_cast<FactorType>(1) );

         // Fetch data from the eight surrounding neighbors
         // possible indices are: c=center cell, h= center cell + h
         const ResType ccc = field_.get( xMin    , yMin    , zMin     , f);
         const ResType hcc = field_.get( xMin + 1, yMin    , zMin     , f);
         const ResType chc = field_.get( xMin    , yMin + 1, zMin     , f);
         const ResType cch = field_.get( xMin    , yMin    , zMin + 1 , f);
         const ResType hhc = field_.get( xMin + 1, yMin + 1, zMin     , f);
         const ResType hch = field_.get( xMin + 1, yMin    , zMin + 1 , f);
         const ResType chh = field_.get( xMin    , yMin + 1, zMin + 1 , f);
         const ResType hhh = field_.get( xMin + 1, yMin + 1, zMin + 1 , f);

         static const FactorType ONE = 1;

         return          xDiff   *         yDiff   *         zDiff   * hhh
               + ( ONE - xDiff ) * ( ONE - yDiff ) * ( ONE - zDiff ) * ccc
               + ( ONE - xDiff ) *         yDiff   *         zDiff   * chh
               +         xDiff   * ( ONE - yDiff ) *         zDiff   * hch
               +         xDiff   *         yDiff   * ( ONE - zDiff ) * hhc
               + ( ONE - xDiff ) * ( ONE - yDiff ) *         zDiff   * cch
               +         xDiff   * ( ONE - yDiff ) * ( ONE - zDiff ) * hcc
               + ( ONE - xDiff ) *         yDiff   * ( ONE - zDiff ) * chc ;
      }

      inline ResType operator()( const Vector3< real_t > & p, cell_idx_t f ) const { return (*this)( p[0], p[1], p[2], f ); }
      inline ResType operator()( const Vector3< real_t > & p )               const { return (*this)( p[0], p[1], p[2] );    }

      const AABB & getValidRegion() const {
         return validRegion_;
      }

   protected:
      const GlField & field_;
      AABB validRegion_;
   };




} // namespace field
} // namespace walberla




