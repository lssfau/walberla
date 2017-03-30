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
//! \file ImageTest.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "geometry/structured/RGBAImage.h"
#include "geometry/structured/GrayScaleImage.h"

#include <iostream>

using namespace walberla;


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );


   using geometry::RGBAImage;

   RGBAImage::pixel_t pixelValue;

   std::string colorStr ("ff1100");
   std::stringstream ss( colorStr );
   ss >> pixelValue;

   WALBERLA_CHECK_EQUAL( pixelValue.color.r, 0xff );
   WALBERLA_CHECK_EQUAL( pixelValue.color.g, 0x11 );
   WALBERLA_CHECK_EQUAL( pixelValue.color.b, 0x00 );
   WALBERLA_CHECK_EQUAL( pixelValue.color.a, 0xff );

   RGBAImage img( 100,100 );
   img.setToWhite();

   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.r , 0xff );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.g , 0xff );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.b , 0xff );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.a , 0xff );

   img.getElement(0,0) = pixelValue;

   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.r , 0xff );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.g , 0x11 );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.b , 0x00 );
   WALBERLA_CHECK_EQUAL( img.getElement(0,0).color.a , 0xff );

   return 0;
}
