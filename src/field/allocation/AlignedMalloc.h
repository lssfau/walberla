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
//! \file AlignedMalloc.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief (De) Allocation of aligned memory
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <cstdlib>


namespace walberla {
namespace field {

   //*******************************************************************************************************************
   /*!
    * Allocates memory starting at an aligned location
    *
    * \ingroup field
    *
    * Uses standard malloc, allocates slightly more memory than needed to
    * ensure alignment. Memory allocated with aligned_malloc can only be freed
    * with aligned_free()
    *
    * \param size      size of allocated memory, total allocation will be greater by
    *                  'alignment' + sizeof(void*) bytes
    * \param alignment  the alignment boundary ,
    *                   **important** has to be a power of 2!   */
   //*******************************************************************************************************************
   void *aligned_malloc( uint_t size, uint_t alignment );



   //*******************************************************************************************************************
   /*!
    * Allocates memory such that (ptr+offset) is aligned
    *
    * \ingroup field
    *
    * This function is useful when allocating aligned memory for data with ghost layers:
    * x(-1) , x(0),  x(1),  x(2),  x(3), ....
    * and x(0) should be aligned. In this case offset should be chosen as nrOfGhostLayers*sizeof( elementType )
    *
    * \param size       see aligned_malloc()
    * \param alignment  see aligned_malloc()
    * \param offset     offset in bytes such that (resulting pointer + offset) is aligned
    * */
   //*******************************************************************************************************************
   void *aligned_malloc_with_offset( uint_t size, uint_t alignment, uint_t offset );


   /****************************************************************************************************************//**
    * Analogous to free for memory allocated with aligned_malloc
    *
    * \ingroup field
    *
    * \param ptr  The pointer returned by aligned_malloc
    *******************************************************************************************************************/
   void aligned_free( void *ptr );


} // namespace field
} // namespace walberla


