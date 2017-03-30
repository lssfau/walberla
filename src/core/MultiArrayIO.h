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
//! \file MultiArrayIO.h
//! \ingroup config
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once


#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning( disable : 4189 )
#  pragma warning( disable : 4100 )
#  pragma warning( disable : 4458 )
#  pragma warning( disable : 4459 )
#  pragma warning( disable : 4510 )
#  pragma warning( disable : 4610 )
#endif //_MSC_VER

#include <boost/multi_array.hpp>

#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER


namespace boost {

   // 1D Arrays
   template< typename T>
   std::istream & operator>> ( std::istream & is, boost::multi_array<T,1> & arr );

   template<typename T>
   std::ostream & operator<< ( std::ostream & os, const boost::multi_array<T,1> & arr );



   // 2D Arrays
   template<typename T>
   std::istream & operator>> ( std::istream & is, boost::multi_array<T,2> & arr );

   template<typename T>
   std::ostream & operator<< ( std::ostream & os, const boost::multi_array<T,2> & arr );



} // namespace walberla


#include "MultiArrayIO.impl.h"
