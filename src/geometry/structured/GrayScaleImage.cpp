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
//! \file GrayScaleImage.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "GrayScaleImage.h"
#include "core/Abort.h"

#include "lodepng.h"

#include <sstream>

namespace walberla   {
namespace geometry   {

   GrayScaleImage::GrayScaleImage( uint_t _width, uint_t _height )
   {
      size_[0] = _width;
      size_[1] = _height;
      image_.resize( size_[0] * size_[1] );
   }

   GrayScaleImage::GrayScaleImage( const std::string & pngFilename )
   {
      unsigned int tmpWidth;
      unsigned int tmpHeight;
      unsigned int error = lodepng::decode( image_, tmpWidth, tmpHeight, pngFilename, LCT_GREY, 8 );
      size_[0] = tmpWidth;
      size_[1] = tmpHeight;

      if ( error )
         WALBERLA_ABORT( "Error while loading PNG file:\n" << lodepng_error_text(error)  );
   }

   void GrayScaleImage::save( const std::string & pngFilename )
   {
      uint32_t error = lodepng::encode( pngFilename, image_,
                                        uint32_c( size_[0] ), uint32_c( size_[1] ),
                                        LCT_GREY, 8 );

      if ( error )
         WALBERLA_ABORT( "Error while loading PNG file:\n" << lodepng_error_text(error)  );
   }


   real_t GrayScaleImage::operator() ( cell_idx_t x, cell_idx_t y ) const
   {
      static const real_t maxVal = real_c( std::numeric_limits<unsigned char>::max() );
      return real_c ( getElement(x, y) ) / maxVal;
   }

   void GrayScaleImage::setElement( cell_idx_t x, cell_idx_t y, real_t val)
   {
      WALBERLA_ASSERT_LESS_EQUAL   ( val, real_t(1.0) );
      WALBERLA_ASSERT_GREATER_EQUAL( val, real_t(0.0) );
      getElement(x,y) = static_cast<unsigned char>( real_c( std::numeric_limits<unsigned char>::max() ) * val  );
   }

   GrayScaleImage GrayScaleImage::getResizedImage( uint_t newWidth, uint_t newHeight, bool bilinear ) const
   {
      if ( newWidth == size_[0]  && newHeight == size_[1] )
         return *this;


      GrayScaleImage resizedImage;

      resizedImage.size_[0] = newWidth;
      resizedImage.size_[1] = newHeight;

      resizedImage.image_.resize( newWidth * newHeight );

      if ( bilinear )
      {
         real_t scaleX = real_c( size_[0]-1 ) / real_c( newWidth );
         real_t scaleY = real_c( size_[1]-1 ) / real_c( newHeight);

         for( cell_idx_t y = 0; y < cell_idx_c( newHeight ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( newWidth ); ++x )
            {
               real_t oldX = real_c(x) * scaleX;
               real_t oldY = real_c(y) * scaleY;
               cell_idx_t oldXi = cell_idx_c( oldX );
               cell_idx_t oldYi = cell_idx_c( oldY );
               real_t xDiff = oldX - real_c(oldXi);
               real_t yDiff = oldY - real_c(oldYi);

               // bilinear interpolation

               resizedImage.getElement( x, y ) =
                        uint8_c(
                        (1 - xDiff) * (1 - yDiff ) * getElement( oldXi    , oldYi    ) +
                             xDiff  * (1 - yDiff ) * getElement( oldXi + 1, oldYi    ) +
                        (1 - xDiff) *      yDiff   * getElement( oldXi    , oldYi + 1) +
                             xDiff  *      yDiff   * getElement( oldXi + 1, oldYi + 1) );
            }
      }
      else
      {
         real_t scaleX = real_c( size_[0] ) / real_c( newWidth );
         real_t scaleY = real_c( size_[1] ) / real_c( newHeight);

         for( cell_idx_t y = 0; y < cell_idx_c( newHeight ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( newWidth ); ++x )
            {
               real_t oldX = real_c(x) * scaleX;
               real_t oldY = real_c(y) * scaleY;
               cell_idx_t oldXi = cell_idx_c( oldX );
               cell_idx_t oldYi = cell_idx_c( oldY );

               resizedImage.getElement( x, y ) = getElement( oldXi, oldYi );
            }
      }

      return resizedImage;
   }

   GrayScaleImage::pixel_t GrayScaleImage::pixelValueFromString( const std::string & str )
   {
      std::stringstream ss (str);
      int tmp;
      ss >> std::hex >> tmp;
      return uint8_c( tmp );
   }

} // namespace geometry
} // namespace walberla









