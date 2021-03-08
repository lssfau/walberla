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
//! \file GrayScaleImage.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <string>
#include <vector>


namespace walberla   {
namespace geometry   {


   //*******************************************************************************************************************
   /*!
   * A gray scale image
   *
   * \ingroup geometry
   *
   * Can be loaded from png files
   */
   //*******************************************************************************************************************
   class GrayScaleImage
   {
   public:
      using pixel_t = unsigned char;

      GrayScaleImage( uint_t _width, uint_t _height );
      GrayScaleImage( const std::string & pngFilename );


      void save( const std::string & pngFilename );

      GrayScaleImage getResizedImage( uint_t newWidth, uint_t newHeight, bool bilinear=true ) const;

      uint_t width()  const { return size_[0]; }
      uint_t height() const { return size_[1]; }

      uint_t size( uint_t coord ) const;

      real_t operator() ( cell_idx_t x, cell_idx_t y ) const;

      void setElement( cell_idx_t x, cell_idx_t y, real_t val);

      unsigned char & getElement ( cell_idx_t x, cell_idx_t y );
      unsigned char   getElement ( cell_idx_t x, cell_idx_t y ) const;


      static pixel_t pixelValueFromString( const std::string & str );

   protected:
      GrayScaleImage() = default;

      uint_t size_[2];                   //< 0=width,  1=height
      std::vector<unsigned char> image_; //< raw pixels
   };



   inline unsigned char &  GrayScaleImage::getElement ( cell_idx_t x, cell_idx_t y ) {
      WALBERLA_ASSERT_GREATER_EQUAL( x, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( y, 0 );
      WALBERLA_ASSERT_LESS( x, cell_idx_c( size_[0] ) );
      WALBERLA_ASSERT_LESS( y, cell_idx_c( size_[1] ) );
      const uint_t yFlip = size_[1] - uint_c(y) - uint_t(1);
      return image_[ yFlip * size_[0] + uint_c(x) ];
   }

   inline unsigned char GrayScaleImage::getElement ( cell_idx_t x, cell_idx_t y ) const {
      return const_cast<GrayScaleImage*> ( this )->getElement(x,y);
   }


   inline uint_t GrayScaleImage::size( uint_t coord ) const
   {
      WALBERLA_ASSERT_LESS( coord, 2);
      return size_[coord];
   }




} // namespace geometry
} // namespace walberla


