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
//! \file MetisWrapper.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h" // for macro WALBERLA_BUILD_WITH_METIS

#ifdef _MSC_VER
   #pragma push_macro( "INT32_MIN" )
   #pragma push_macro( "INT32_MAX" )
   #pragma push_macro( "INT64_MIN" )
   #pragma push_macro( "INT64_MAX" )

   #ifdef INT32_MIN
      #undef INT32_MIN
   #endif

   #ifdef INT32_MAX
      #undef INT32_MAX
   #endif

   #ifdef INT64_MIN
      #undef INT64_MIN
   #endif

   #ifdef INT64_MAX
      #undef INT64_MAX
   #endif
#endif

#ifdef WALBERLA_BUILD_WITH_METIS
// external software includes
#include "metis.h"
#endif

#ifdef _MSC_VER
#pragma pop_macro( "INT64_MAX" )
#pragma pop_macro( "INT64_MIN" )
#pragma pop_macro( "INT32_MAX" )
#pragma pop_macro( "INT32_MIN" )
#endif
