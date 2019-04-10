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
//! \file BoundaryFromDomainBorder.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundarySetter.h"

#include "Initializer.h"

#include "boundary/Boundary.h"

#include "core/Abort.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "stencil/D3Q7.h"
#include "stencil/Directions.h"


namespace walberla {
namespace geometry {
namespace initializer {



//**********************************************************************************************************************
/*! Sets a boundary condition at a certain border of the domain (borders are specified by directions).
*
* Configuration file syntax:
   \verbatim
     <initializerID>
     {
         direction    W;         // can be W,E,S,N,T,B  and all
         walldistance 0;         // optional (default to 0),
                                 // 0 means the first interior layer is set, -1 means ghost layer
                                 // 1 first inner slice ...

        ghostLayersToInitialize  // how many ghost layers are initialized in the direction normal to the border
                                 // Example: when initializing northern boundary, the parameter controls how far
                                 // the slice goes into the W,E,T,B directions, i.e. how many ghost layers are
                                 // initialized in these directions
                                 // If this parameter is greater than the existing number of ghost layers,
                                 // no error is reported, but all ghostlayers are initialized.
                                 // Do not mix up this parameter with the parameter "walldistance" which
                                 // specifies the distance of the slice to the (in this example) northern boundary.

         <FlagUID>               {  <boundaryConfiguration>   }
     }
   \endverbatim

*  The order of these blocks in the configuration file is important, when regions overlap.
*  later config blocks overwrite the previous ones.
*
*/
//**********************************************************************************************************************
template<typename BoundaryHandling_T>
class BoundaryFromDomainBorder : public Initializer
{
public:
   BoundaryFromDomainBorder( StructuredBlockStorage & blocks, BlockDataID handlerBlockDataID,
                             CellInterval globalDomain = CellInterval() );

   void init ( const Config::BlockHandle & blockHandle );

   virtual void init( BlockStorage & blockStorage, const Config::BlockHandle & blockHandle );


   void init( FlagUID flagUID, stencil::Direction direction,
              cell_idx_t wallDistance = 0, cell_idx_t ghostLayersToInitialize = std::numeric_limits<cell_idx_t>::max() );

   void init( BoundaryUID boundaryUID, stencil::Direction direction, const shared_ptr<BoundaryConfiguration> & conf,
              cell_idx_t wallDistance = 0,
              cell_idx_t ghostLayersToInitialize = std::numeric_limits<cell_idx_t>::max() );

   void initAllBorders( FlagUID flagUID, cell_idx_t wallDistance = 0,
                        cell_idx_t ghostLayersToInitialize = std::numeric_limits<cell_idx_t>::max() );

   void initAllBorders( BoundaryUID boundaryUID, const shared_ptr<BoundaryConfiguration> & conf,
                        cell_idx_t wallDistance = 0, cell_idx_t ghostLayersToInitialize = std::numeric_limits<cell_idx_t>::max() );


protected:

   void init( stencil::Direction direction,
              cell_idx_t wallDistance,
              cell_idx_t ghostLayersToInitialize,
              BoundarySetter<BoundaryHandling_T> & boundarySetter );

   StructuredBlockStorage & blocks_;
   BlockDataID handlerBlockDataID_;

   CellInterval globalDomain_;
};




} // namespace initializer
} // namespace geometry
} // namespace walberla


#include "BoundaryFromDomainBorder.impl.h"

