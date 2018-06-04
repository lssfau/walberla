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
//! \file MakeBlockDataInitFunction.h
//! \ingroup domain_decomposition
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "IBlock.h"

#include <functional>


namespace walberla {
namespace domain_decomposition {



/// \cond internal
namespace internal
{
   template<class T, typename... Args>
   T * newFunc( const IBlock* const, Args&... args ) {
      return new T(args...);
   }
} // namespace internal
/// \endcond



template<class T, typename... Args>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(Args&&... args) {
   return std::bind( internal::newFunc<T,Args...>, std::placeholders::_1, std::forward<Args>(args)... );
}



} // namespace domain_decomposition

using domain_decomposition::makeBlockDataInitFunction;

} // namespace walberla
