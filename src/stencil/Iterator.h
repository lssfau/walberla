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
//! \file Iterator.h
//! \ingroup stencil
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Directions.h"
#include "core/debug/Debug.h"

namespace walberla {
namespace stencil {


//**********************************************************************************************************************
/*! \brief Iterator over all directions contained in a stencil
*
* \ingroup stencil
*
* See Stencil.in.h for documentation how to use the iterator
*
************************************************************************************************************************/
template<typename Stencil>
class Iterator
{
public:
   explicit Iterator(uint_t i) : i_(i) {}

   //** Operators  *****************************************************************************************************
   /*! \name Operators*/
   //@{
   inline Iterator & operator++()                   { ++i_;  return *this; }
   inline bool operator==(const Iterator & o) const { return i_ == o.i_; }
   inline bool operator!=(const Iterator & o) const { return i_ != o.i_; }
   //}
   //*******************************************************************************************************************


   //** Access Functions   *********************************************************************************************
   /*! \name Access Functions*/
   //@{
   inline Direction          operator*()         const { return Stencil::dir[i_];                       }
   inline Direction          direction()         const { return Stencil::dir[i_];                       }
   inline Direction          inverseDir()        const { return stencil::inverseDir[Stencil::dir[i_]];  }
   inline uint_t             toIdx()             const { return i_;                                     }
   inline uint_t             toInvIdx()          const { return Stencil::idx[inverseDir()];             }
   inline int                cx()                const { return stencil::cx[Stencil::dir[i_]];          }
   inline int                cy()                const { return stencil::cy[Stencil::dir[i_]];          }
   inline int                cz()                const { return stencil::cz[Stencil::dir[i_]];          }
   inline int                c( const uint_t d ) const { WALBERLA_ASSERT_LESS( d, 3 );
                                                         return stencil::c[d][Stencil::dir[i_]];        }
   inline real_t             length()            const { return stencil::dirLength[Stencil::dir[i_]];   }
   inline BinaryDirection    binaryDir()         const { return stencil::dirToBinary[Stencil::dir[i_]]; }
   inline const std::string& dirString()         const { return stencil::dirToString[Stencil::dir[i_]]; }
   inline Direction          mirrorX()           const { return stencil::mirrorX[Stencil::dir[i_]];     }
   inline Direction          mirrorY()           const { return stencil::mirrorY[Stencil::dir[i_]];     }
   inline Direction          mirrorZ()           const { return stencil::mirrorZ[Stencil::dir[i_]];     }
   //}
   //*******************************************************************************************************************

private:
   uint_t i_;
};


} // namespace stencil
} // namespace walberla
