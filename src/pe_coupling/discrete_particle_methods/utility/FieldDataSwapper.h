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
//! \file FieldDataSwapper.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Swaps the data of two identical fields.
 *
 * This functionality is e.g. used to store a former version of the velocity field to calculate the time derivative
 * of the velocity numerically.
 */
template< typename Field_T >
class FieldDataSwapper
{
public:

   FieldDataSwapper( const BlockDataID & srcFieldID, const BlockDataID & dstFieldID )
      : srcFieldID_( srcFieldID ), dstFieldID_( dstFieldID )
   { }

   void operator()(IBlock * const block)
   {
      Field_T* srcField = block->getData< Field_T >( srcFieldID_ );
      Field_T* dstField = block->getData< Field_T >( dstFieldID_ );
      srcField->swapDataPointers(dstField);
   }

private:
   const BlockDataID srcFieldID_;
   const BlockDataID dstFieldID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
