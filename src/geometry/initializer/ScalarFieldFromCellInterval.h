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
//! \file ScalarFieldFromCellInterval.h
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "geometry/initializer/Initializer.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "core/math/Parser.h"


namespace walberla {
namespace geometry {
namespace initializer {


//*******************************************************************************************************************
/*!
* Initializes a scalar field using a CellInterval
*
* \ingroup geometry
*
* Example:
  \verbatim
     <InitializerId> {
        CellInterval [ (3,1,5) ... (4,5,8) ];
        id           0;    // If given a vector of scalar fields, which one of them to operate on.
                           // Should be zero (the default value) if given a scalar field directly.
        value        3.14; // value to set on each cell in the interval
     }
  \endverbatim

or
  \verbatim
     <InitializerId> {
        min < 3,1,5>;
        max < 4,5,8>;
        id           0;
        value        3.14;
     }
  \endverbatim
*
*/
//*******************************************************************************************************************
template <typename Field_T>
class ScalarFieldFromCellInterval : public Initializer
{
public:

   typedef typename Field_T::value_type Value_T;

   ScalarFieldFromCellInterval( StructuredBlockStorage & blocks, BlockDataID fieldId )
      : ScalarFieldFromCellInterval(blocks, std::vector<BlockDataID>(1, fieldId))
   {}
   
   ScalarFieldFromCellInterval( StructuredBlockStorage & blocks, std::vector<BlockDataID> fieldId );



   /*************************************************************************************************************//**
   * Initializes the scalar field using parameters of config block
   *****************************************************************************************************************/
   virtual void init ( BlockStorage & blockStorage, const Config::BlockHandle & block );
   void init( const Config::BlockHandle & blockHandle );


   /*************************************************************************************************************//**
   * Function for manually setting a scalar field on a CellInterval
   *
   * \param interval             the cell interval
   * \param value                which value to set in the field for all cells in the interval
   * \param id                   which field to operate on (if operating on a vector of fields), defaults to 0
   *****************************************************************************************************************/
   void init( const CellInterval & interval, Value_T value, std::vector<BlockDataID>::size_type id = 0 );
   /*************************************************************************************************************//**
   * Function for manually setting a scalar field on a CellInterval
   *
   * \param interval             the cell interval
   * \param parser               a function parser which will have the variables x,y,z bound before it is evaluated
   * \param id                   which field to operate on (if operating on a vector of fields), defaults to 0
   *****************************************************************************************************************/
   void init( const CellInterval & interval, math::FunctionParser & parser, std::vector<BlockDataID>::size_type id = 0 );


protected:

   StructuredBlockStorage & structuredBlockStorage_;
   std::vector<BlockDataID> scalarFieldID_;

};


} // namespace initializer
} // namespace geometry
} // namespace walberla

#include "ScalarFieldFromCellInterval.impl.h"
