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
//! \file BoundaryFromImage.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Initializes boundary handling using an image
//
//======================================================================================================================

#pragma once

#include "BoundarySetter.h"
#include "Initializer.h"

#include "boundary/Boundary.h"

#include "core/Abort.h"
#include "core/config/Config.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagUID.h"

#include <map>

namespace walberla {
namespace geometry {
namespace initializer {

   //*******************************************************************************************************************
   /*!
   * Sets boundary conditions according to a gray scale image.
   *
   * \ingroup geometry
   *
   *
   *
     Example:
     \verbatim
        <InitializerId> {
           file                 domain.png;
           extrusionCoordinate  1;  // defaults to 2
           lowerExtrusionLimit  5   // defaults to 0
           upperExtrusionLimit  10; // defaults to maximum of extrusionCoordinate

           xOffset              5; // "x" offset of the image in the domain
           yOffset              6; // "y" offset of the image in the domain
           // or
           rescaleToDomain  true;

           BoundaryValueMapping {
              <BoundaryUID> { <boundary Config> } // boundary name i.e. name of boundary flag
              #or
              flag <FlagUID>;

              value      ff0000;                  // boundary is set in cells where image has exactly
                                                  // this value ( here: red)
                                                  // value has to be specified in hex notation
                                                  // for rgba images the format is ff11aaff
                                                  // ( red=255, green=11, blue=aa, alpha=ff)
                                                  // where the alpha channel is optional
                                                  // for grayscale images only the gray value has to be specified:
                                                  // for example "ff" for white
           }

           BoundaryValueMapping {
              NoSlip     {}
              value      125;
           }
        }

     \endverbatim
   */
   //*******************************************************************************************************************
   template<typename BoundaryHandling_T, typename Image_T >
   class BoundaryFromImage : public Initializer
   {
   public:
      typedef typename Image_T::pixel_t pixel_t;
      typedef std::map<pixel_t, BoundarySetter<BoundaryHandling_T> > BoundarySetters;


      BoundaryFromImage( StructuredBlockStorage & blocks, BlockDataID handlerBlockDataID );

      /*************************************************************************************************************//**
      * Initializes the scalar field using parameters of config block
      * for syntax see class documentation
      *****************************************************************************************************************/
      void init ( BlockStorage & , const Config::BlockHandle & block ) override { init(block); }
      void init( const Config::BlockHandle & block );


      /*************************************************************************************************************//**
      * Initializes boundary values using an image
      *
      * \param img                  the image
      * \param extrusionCoord       the image is set in slices where extrusionCoord is constant
      * \param boundarySetters      a map from pixel value to boundary information
      * \param xOffset              offset for the x IMAGE coordinate ( not cell coordinates )
      * \param yOffset              offset for the y IMAGE coordinate ( not cell coordinates )
      * \param lowerExtrusionLimit  only cells where extrusionCoord is bigger than this parameter are modified
      * \param upperExtrusionLimit  only cells where extrusionCoord is smaller than this parameter are modified
      *****************************************************************************************************************/
      void init( const Image_T & img, uint_t extrusionCoord,
                 BoundarySetters & boundarySetters,
                 cell_idx_t xOffset, cell_idx_t yOffset,
                 cell_idx_t lowerExtrusionLimit,
                 cell_idx_t upperExtrusionLimit );


   protected:
      Image_T resizeImageToDomain( const Image_T & img, uint_t extrusionCoord ) const;


      StructuredBlockStorage & structuredBlockStorage_;
      BlockDataID              handlingID_;
   };



} // namespace initializer
} // namespace geometry
} // namespace walberla

#include "BoundaryFromImage.impl.h"
