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
   template<class T>
   T * newFunc( const IBlock* const ) {
      return new T();
   }

   template<class T, class P1>
   T * newFunc( const IBlock* const , const P1 & p1 ) {
      return new T(p1);
   }

   template<class T, class P1, class P2>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2 ) {
      return new T(p1, p2);
   }

   template<class T, class P1, class P2, class P3>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3 ) {
      return new T(p1, p2, p3);
   }

   template<class T, class P1, class P2, class P3, class P4>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4 ) {
      return new T(p1, p2, p3, p4);
   }

   template<class T, class P1, class P2, class P3, class P4, class P5>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5 ) {
      return new T(p1, p2, p3, p4, p5);
   }

   template<class T, class P1, class P2, class P3, class P4, class P5, class P6>
   T * newFunc( const IBlock* const ,  const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6 ) {
      return new T(p1, p2, p3, p4, p5, p6);
   }

   template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7 ) {
      return new T(p1, p2, p3, p4, p5, p6,p7);
   }

   template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7, class P8>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7, const P8 & p8 ) {
      return new T(p1, p2, p3, p4, p5, p6, p7,p8 );
   }

   template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7, class P8, class P9>
   T * newFunc( const IBlock* const , const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7, const P8 & p8, const P9 & p9 ) {
      return new T(p1, p2, p3, p4, p5, p6, p7,p8, p9 );
   }

} // namespace internal
/// \endcond



template<class T>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction() {
   return std::bind( internal::newFunc<T>,  std::placeholders::_1);
}

template<class T, class P1>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1) {
   return std::bind( internal::newFunc<T,P1>,  std::placeholders::_1,p1);
}

template<class T, class P1, class P2>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2) {
   return std::bind( internal::newFunc<T,P1,P2>,  std::placeholders::_1,p1,p2);
}

template<class T, class P1, class P2, class P3>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3) {
   return std::bind( internal::newFunc<T,P1,P2,P3>,  std::placeholders::_1,p1,p2,p3);
}

template<class T, class P1, class P2, class P3, class P4>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4>,  std::placeholders::_1,p1,p2,p3,p4);
}

template<class T, class P1, class P2, class P3, class P4, class P5>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4,P5>,  std::placeholders::_1,p1,p2,p3,p4,p5);
}

template<class T, class P1, class P2, class P3, class P4, class P5, class P6>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4,P5,P6>,  std::placeholders::_1,p1,p2,p3,p4,p5,p6);
}

template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4,P5,P6,P7>,  std::placeholders::_1,p1,p2,p3,p4,p5,p6,p7);
}

template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7, class P8>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7, const P8 & p8) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4,P5,P6,P7,P8>,  std::placeholders::_1,p1,p2,p3,p4,p5,p6,p7,p8);
}

template<class T, class P1, class P2, class P3, class P4, class P5, class P6, class P7, class P8, class P9>
std::function< T* ( const IBlock* const block ) >
makeBlockDataInitFunction(const P1 & p1, const P2 & p2, const P3 & p3, const P4 & p4, const P5 & p5, const P6 & p6, const P7 & p7, const P8 & p8, const P9 & p9) {
   return std::bind( internal::newFunc<T,P1,P2,P3,P4,P5,P6,P7,P8,P9>,  std::placeholders::_1,p1,p2,p3,p4,p5,p6,p7,p8,p9);
}



} // namespace domain_decomposition

using domain_decomposition::makeBlockDataInitFunction;

} // namespace walberla
