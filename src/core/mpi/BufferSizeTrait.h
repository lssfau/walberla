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
//! \file BufferSizeTrait.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include "core/DataTypes.h"

#include <type_traits>


namespace walberla {
namespace mpi {

#ifdef WALBERLA_BUFFER_DEBUG
   const uint_t BUFFER_DEBUG_OVERHEAD = 2;
#else
   const uint_t BUFFER_DEBUG_OVERHEAD = 0;
#endif

   template< typename T, class Enable = void >
   struct BufferSizeTrait   { };

   #define BufferSizeTraitSpecialization( TYPE ) \
   template<>\
   struct BufferSizeTrait<TYPE>\
   {\
      static const bool constantSize = true;\
      static const size_t size = sizeof( TYPE ) + BUFFER_DEBUG_OVERHEAD;\
   };

   BufferSizeTraitSpecialization(bool)
   BufferSizeTraitSpecialization(char)
   BufferSizeTraitSpecialization(short)
   BufferSizeTraitSpecialization(int)
   BufferSizeTraitSpecialization(long)
   BufferSizeTraitSpecialization(long long)
   BufferSizeTraitSpecialization(signed char)
   BufferSizeTraitSpecialization(unsigned char)
   BufferSizeTraitSpecialization(unsigned short)
   BufferSizeTraitSpecialization(unsigned int)
   BufferSizeTraitSpecialization(unsigned long)
   BufferSizeTraitSpecialization(unsigned long long)
   BufferSizeTraitSpecialization(float)
   BufferSizeTraitSpecialization(double)

   #undef BufferSizeTraitSpecialization

   // Specialization for all enums
   template< typename T>
   struct BufferSizeTrait< T, typename std::enable_if< std::is_enum< T >::value >::type >
   {
      static const bool constantSize = true;
      static const size_t size = sizeof( T ) + BUFFER_DEBUG_OVERHEAD;
   };

} // namespace mpi
} // namespace walberla


