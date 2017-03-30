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
//! \file FlagFunctions.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Bit operations on integer data types.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace field {

//** Operations on bit-masks    ****************************************************************************************
/*! \name Operations on bit-masks
 *  \ingroup field  */
//@{
template<class T> inline void addMask         (      T & val, T mask) { static_assert_int_t<T>();          val  = T( val | mask ); }
template<class T> inline void removeMask      (      T & val, T mask) { static_assert_int_t<T>();          val &= T(~mask);        }
template<class T> inline bool isMaskSet       (const T & val, T mask) { static_assert_int_t<T>(); return ( val & mask) == mask;    }
template<class T> inline bool isPartOfMaskSet (const T & val, T mask) { static_assert_int_t<T>(); return ( val & mask) > T(0);     }
//@}

//**********************************************************************************************************************



//** Operations on bit-masks using flags *******************************************************************************
/*! \name Operations on bit-masks using flags
 * \ingroup field
 *
 * Flag is defined as a mask, where exactly one bit is set. This restriction is
 * enforced in debug mode using asserts. In release mode the following functions are equivalent
 * to their *Mask* variants.
 */
//@{
#define F_MSG "Second parameter f=" << f << " has to be a flag (only one bit is 1)"
template<class T> inline bool isFlag          (             T f) { static_assert_int_t<T>(); return f && !( f & (f - T(1)) );                       }
template<class T> inline void addFlag         (      T & v, T f) { static_assert_int_t<T>(); WALBERLA_ASSERT( isFlag(f), F_MSG ); addMask(v,f);          }
template<class T> inline void removeFlag      (      T & v, T f) { static_assert_int_t<T>(); WALBERLA_ASSERT( isFlag(f), F_MSG ); removeMask(v,f);       }
template<class T> inline bool isFlagSet       (const T & v, T f) { static_assert_int_t<T>(); WALBERLA_ASSERT( isFlag(f), F_MSG ); return isMaskSet(v,f); }
#undef F_MSG
//@}


} // namespace field
} // namespace walberla





//======================================================================================================================
//
//  EXPORTS
//
//======================================================================================================================

namespace walberla {
   using field::addMask;
   using field::removeMask;
   using field::isMaskSet;
   using field::isPartOfMaskSet;

   using field::addFlag;
   using field::removeFlag;
   using field::isFlagSet;
}
