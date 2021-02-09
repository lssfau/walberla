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
//! \file BoundaryFromCellInterval.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundarySetter.h"
#include "Initializer.h"

#include "boundary/Boundary.h"
#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "domain_decomposition/BlockStorage.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <string>


namespace walberla {
namespace geometry {
namespace initializer {



//*******************************************************************************************************************
/*! Initializes a boundary handler using a CellInterval
*   (works for both boundary handlings and boundary handling collections).
*
* Configuration file syntax:
 \verbatim
     <initializerID>
     {
         CellInterval [ (3,1,5) ... (4,5,8) ]
         <FlagUID> {  <boundaryConfiguration>   }
     }
 \endverbatim

 or
 \verbatim
     <initializerID>
     {
        min < 3,1,5>;
        max < 4,5,8>;
        <FlagUID> {  <boundaryConfiguration>   }
     }
 \endverbatim

*
* \ingroup geometry
*/
//*******************************************************************************************************************
template <typename BoundaryHandlerT>
class BoundaryFromCellInterval : public Initializer
{
public:
   BoundaryFromCellInterval( StructuredBlockStorage & structuredBlockStorage, BlockDataID & boundaryHandlerID );


   /*************************************************************************************************************//**
   * Sets the boundary using information from a config-block
   *
   * The cell interval and the parameters for the boundary are read from configuration.
   *
   *****************************************************************************************************************/
   void init( BlockStorage & , const Config::BlockHandle & blockHandle ) override { return init(blockHandle); }
   void init( const Config::BlockHandle & blockHandle );


   /*************************************************************************************************************//**
   * Function for manually setting a boundary condition on a CellInterval
   *****************************************************************************************************************/

   void init( const CellInterval & ci, const BoundaryUID & uid, const shared_ptr<BoundaryConfiguration> & bcConfig );
   void init( const CellInterval & ci, const FlagUID & uid );


protected:
   void init( const CellInterval & globalCellInterval, BoundarySetter<BoundaryHandlerT> & boundarySetter );

   StructuredBlockStorage & structuredBlockStorage_;
   BlockDataID              boundaryHandlerID_;
};



} // namespace initializer
} // namespace geometry
} // namespace walberla

#include "BoundaryFromCellInterval.impl.h"
