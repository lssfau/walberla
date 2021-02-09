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
//! \file ScalarFieldFromGrayScaleImage.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Initializer.h"
#include "geometry/structured/GrayScaleImage.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace geometry {
namespace initializer {


   //*******************************************************************************************************************
   /*!
   * Initializes a scalar field, using the pixel values from a gray scale image
   * The `real_t` field is initialized with values between 0 and 1
   *
   * \ingroup geometry
   *
     Example:
     \verbatim
        <InitializerId> {
           file                 domain.png;
           extrusionCoordinate  1; // y coordinate
           lowerExtrusionLimit  5;
           upperExtrusionLimit  10;
           xOffset              5;
           yOffset              6;
        }
     \endverbatim
   *
   */
   //*******************************************************************************************************************
   class ScalarFieldFromGrayScaleImage : public Initializer
   {
   public:

      ScalarFieldFromGrayScaleImage( StructuredBlockStorage & blocks, BlockDataID scalarFieldID );



      /*************************************************************************************************************//**
      * Initializes the scalar field using parameters of config block
      * for syntax see class documentation
      *****************************************************************************************************************/
      void init ( BlockStorage & blockStorage, const Config::BlockHandle & block ) override;


      /*************************************************************************************************************//**
      * Initializes scalar field using a gray scale image
      *
      * \param img                  the image
      * \param extrusionCoord       the image is set in slices where extrusionCoord is constant
      * \param xOffset              offset for the x IMAGE coordinate ( not cell coordinates )
      * \param yOffset              offset for the y IMAGE coordinate ( not cell coordinates )
      * \param lowerExtrusionLimit  only cells where extrusionCoord is bigger than this parameter are modified
      * \param upperExtrusionLimit  only cells where extrusionCoord is smaller than this parameter are modified
      *****************************************************************************************************************/
      void init( const GrayScaleImage & img, uint_t extrusionCoord,
                 cell_idx_t xOffset, cell_idx_t yOffset,
                 cell_idx_t lowerExtrusionLimit = std::numeric_limits<cell_idx_t>::min(),
                 cell_idx_t upperExtrusionLimit = std::numeric_limits<cell_idx_t>::max() );



      /*************************************************************************************************************//**
      *  Initializes scalar field using a gray scale image
      *
      *  \param rescaleToDomain   if true  the image is first rescaled to match the size of the domain
      *  other parameters as for function above
      *****************************************************************************************************************/
      void init( const GrayScaleImage & img, uint_t extrusionCoord, bool rescaleToDomain = true,
                 cell_idx_t lowerExtrusionLimit = std::numeric_limits<cell_idx_t>::min(),
                 cell_idx_t upperExtrusionLimit = std::numeric_limits<cell_idx_t>::max() );


   protected:

      StructuredBlockStorage & structuredBlockStorage_;
      BlockDataID              scalarFieldID_;

   };


} // namespace initializer
} // namespace geometry
} // namespace walberla


