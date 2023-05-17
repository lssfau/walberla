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
//! \file AlignedAllocation.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#pragma once

#include "core/DataTypes.h"


namespace walberla {
namespace gpu
{

   void *allocate_aligned_with_offset( uint_t size, uint_t alignment, uint_t offset );


   void free_aligned_with_offset( void *ptr );


   void *allocate_pitched_with_offset( size_t &pitchOut, size_t width, size_t height,
                                       size_t alignment, size_t alignmentOffset );

} // namespace gpu
} // namespace walberla
