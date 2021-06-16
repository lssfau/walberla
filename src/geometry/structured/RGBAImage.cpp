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
//! \file RGBAImage.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "RGBAImage.h"

#include "core/Abort.h"

#include "lodepng.h"

#include <algorithm>
#include <sstream>
#include <iomanip>


namespace walberla   {
namespace geometry   {


   bool RGBAImage::pixel_t::operator< ( const pixel_t & o ) const {
      return std::lexicographical_compare(   values,   values+4, o.values, o.values+4 );
   }

   bool RGBAImage::pixel_t::operator== ( const pixel_t & o ) const {
      return std::equal(   values,   values+4, o.values );
   }


   RGBAImage::RGBAImage( uint_t _width, uint_t _height )
   {
      size_[0] = _width;
      size_[1] = _height;
      image_.resize( 4 * size_[0] * size_[1], static_cast<unsigned char>(255) );
   }

   RGBAImage::RGBAImage( const std::string & pngFilename )
   {
      unsigned int tmpWidth;
      unsigned int tmpHeight;
      unsigned int error = lodepng::decode( image_, tmpWidth, tmpHeight, pngFilename, LCT_RGBA, 8 );
      size_[0] = tmpWidth;
      size_[1] = tmpHeight;

      if ( error )
         WALBERLA_ABORT( "Error while loading PNG file:\n" << lodepng_error_text(error)  );
   }

   void RGBAImage::save( const std::string & pngFilename )
   {
      uint32_t error = lodepng::encode( pngFilename, image_,
                                        uint32_c( size_[0] ), uint32_c( size_[1] ),
                                        LCT_RGBA, 8 );

      if ( error )
         WALBERLA_ABORT( "Error while loading PNG file:\n" << lodepng_error_text(error)  );
   }


   real_t RGBAImage::operator() ( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c ) const
   {
      static const real_t maxVal = real_c( std::numeric_limits<unsigned char>::max() );

      const cell_idx_t yFlip = cell_idx_c( size_[1] ) - y - 1;

      return real_c ( getElement(x, yFlip, c) ) / maxVal;
   }

   void RGBAImage::setElement( cell_idx_t x, cell_idx_t y, RGBAImage::Channel c, real_t val)
   {
      WALBERLA_ASSERT_LESS_EQUAL   ( val, real_t(1.0) );
      WALBERLA_ASSERT_GREATER_EQUAL( val, real_t(0.0) );
      getElement(x,y,c) = static_cast<unsigned char> ( real_c( std::numeric_limits<unsigned char>::max() ) * val );
   }

   RGBAImage RGBAImage::getResizedImage( uint_t newWidth, uint_t newHeight, bool bilinear ) const
   {
      if ( newWidth == size_[0]  && newHeight == size_[1] )
         return *this;


      RGBAImage resizedImage;

      resizedImage.size_[0] = newWidth;
      resizedImage.size_[1] = newHeight;

      resizedImage.image_.resize( 4 * newWidth * newHeight );


      if ( bilinear )
      {
         real_t scaleX = real_c( size_[0] - 1 ) / real_c( newWidth );
         real_t scaleY = real_c( size_[1] - 1 ) / real_c( newHeight);

         for( cell_idx_t y = 0; y < cell_idx_c( newHeight ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( newWidth ); ++x )
               for( int c = 0; c < 4; ++c )
               {
                  real_t oldX = real_c(x) * scaleX;
                  real_t oldY = real_c(y) * scaleY;
                  cell_idx_t oldXi = cell_idx_c( oldX );
                  cell_idx_t oldYi = cell_idx_c( oldY );
                  real_t xDiff = oldX - real_c(oldXi);
                  real_t yDiff = oldY - real_c(oldYi);

                  // bilinear interpolation

                  resizedImage.getElement( x, y, Channel(c) ) =
                           uint8_c(
                           (1 - xDiff) * (1 - yDiff ) * getElement( oldXi    , oldYi    , Channel(c)) +
                                xDiff  * (1 - yDiff ) * getElement( oldXi + 1, oldYi    , Channel(c)) +
                           (1 - xDiff) *      yDiff   * getElement( oldXi    , oldYi + 1, Channel(c)) +
                                xDiff  *      yDiff   * getElement( oldXi + 1, oldYi + 1, Channel(c)) );
               }
      }
      else
      {
         real_t scaleX = real_c( size_[0] ) / real_c( newWidth );
         real_t scaleY = real_c( size_[1] ) / real_c( newHeight);

         for( cell_idx_t y = 0; y < cell_idx_c( newHeight ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( newWidth ); ++x )
               for( int c = 0; c < 4; ++c )
               {
                  real_t oldX = real_c(x) * scaleX;
                  real_t oldY = real_c(y) * scaleY;
                  cell_idx_t oldXi = cell_idx_c( oldX );
                  cell_idx_t oldYi = cell_idx_c( oldY );
                  resizedImage.getElement( x, y, Channel(c) ) = getElement( oldXi, oldYi, Channel(c) );
               }
      }

      return resizedImage;
   }


   RGBAImage::pixel_t RGBAImage::pixelValueFromString( const std::string & value )
   {
      RGBAImage::pixel_t res;

      for( uint_t i=0; i < 3; ++i )
         res.values[i] = 0;
      res.values[3] = 255;

      for( uint_t i=0; i < 4 && i < value.size() / 2; ++i)
      {
         std::string hexPair = value.substr( i*2, 2 );
         std::stringstream ss ( hexPair );
         int tmp;
         ss >> std::hex >> tmp;
         res.values[i] = uint8_c( tmp );
      }
      return res;
   }


   std::ostream & operator<< ( std::ostream & os, const RGBAImage::pixel_t & pixel )
   {
      os << std::hex << std::setw(2) << std::setfill('0') << int_c( pixel.values[0] );
      os << std::hex << std::setw(2) << std::setfill('0') << int_c( pixel.values[1] );
      os << std::hex << std::setw(2) << std::setfill('0') << int_c( pixel.values[2] );
      os << std::hex << std::setw(2) << std::setfill('0') << int_c( pixel.values[3] );
      return os;
   }

   std::istream & operator>> ( std::istream & is, RGBAImage::pixel_t & pixel )
   {
      std::string value;
      is >> value;
      pixel = RGBAImage::pixelValueFromString( value );
      return is;
   }

} // namespace geometry
} // namespace walberla









