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
//! \file RGBAImage.h
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
   class RGBAImage
   {
   public:
      enum Channel {
         R = 0, //< red
         G = 1, //< green
         B = 2, //< blue
         A = 3  //< alpha
      };

      struct RGBAColor
      {
         unsigned char r;
         unsigned char g;
         unsigned char b;
         unsigned char a;
      };
      
      struct pixel_t
      {
         union
         {
            RGBAColor color;
            unsigned char values[4];
            uint32_t value;
         };

         bool operator< ( const pixel_t & o ) const;
         bool operator== ( const pixel_t & o ) const;
      };

      RGBAImage( uint_t _width, uint_t _height );
      RGBAImage( const std::string & pngFilename );

      void save( const std::string & pngFilename );

      RGBAImage getResizedImage( uint_t newWidth, uint_t newHeight, bool bilinear=true ) const;

      uint_t width()  const { return size_[0]; }
      uint_t height() const { return size_[1]; }

      uint_t size( uint_t coord ) const;

      void setToWhite() { std::fill( image_.begin(), image_.end(), static_cast<unsigned char>(255) ); }

      real_t operator() ( cell_idx_t x, cell_idx_t y, Channel c ) const;

      void setElement( cell_idx_t x, cell_idx_t y, Channel c, real_t val );

            unsigned char & getElement ( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c );
      const unsigned char & getElement ( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c ) const;


            pixel_t & getElement( cell_idx_t x, cell_idx_t y );
      const pixel_t & getElement( cell_idx_t x, cell_idx_t y ) const;

      static pixel_t pixelValueFromString( const std::string & str );

   protected:
      RGBAImage() = default;

      uint_t size_[2];                   //< 0=width,  1=height
      std::vector<unsigned char> image_; //< raw pixels
   };

   std::ostream & operator<< ( std::ostream & os, const RGBAImage::pixel_t & pixel );
   std::istream & operator>> ( std::istream & is, RGBAImage::pixel_t & pixel );




   inline unsigned char &  RGBAImage::getElement ( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c ) {
      WALBERLA_ASSERT_GREATER_EQUAL( x, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( y, 0 );
      WALBERLA_ASSERT_LESS( x, cell_idx_c( size_[0] ) );
      WALBERLA_ASSERT_LESS( y, cell_idx_c( size_[1] ) );
      const uint_t yFlip = size_[1] - uint_c(y) - uint_t(1);
      return image_[ uint_t(4) * (yFlip * size_[0] + uint_c(x)) + c ];
   }

   inline const unsigned char & RGBAImage::getElement ( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c ) const {
      return const_cast<RGBAImage*> ( this )->getElement(x,y,c);
   }

   inline RGBAImage::pixel_t & RGBAImage::getElement( cell_idx_t x, cell_idx_t y )
   {
      unsigned char & begin = getElement(x,y,R);
      RGBAImage::pixel_t * pixelPtr = reinterpret_cast<pixel_t*>( &begin );
      return *pixelPtr;
   }

   inline const RGBAImage::pixel_t & RGBAImage::getElement( cell_idx_t x, cell_idx_t y ) const
   {
      const unsigned char & begin = getElement(x,y,R);
      const pixel_t * pixelPtr = reinterpret_cast< const pixel_t*>( &begin );
      return *pixelPtr;
   }

   inline uint_t RGBAImage::size( uint_t coord ) const
   {
      WALBERLA_ASSERT_LESS( coord, 2);
      return size_[coord];
   }


} // namespace geometry
} // namespace walberla


