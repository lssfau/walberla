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
//! \file CellInterval.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "CellInterval.h"
#include "CellSet.h"
#include "CellVector.h"


namespace walberla {
namespace cell {



bool CellInterval::overlaps( const CellSet& cellSet ) const {

   if( empty() )
      return false;

   for( CellSet::const_iterator cell = cellSet.begin(); cell != cellSet.end(); ++cell ) {
      if( this->contains( *cell ) )
         return true;
   }

   return false;
}



bool CellInterval::overlaps( const CellVector& cellVector ) const
{
   if( empty() )
      return false;

   return std::any_of(cellVector.begin(),
                      cellVector.end(),
                      [&](const Cell & cell)
                      {return contains( cell );});
}



} // namesapce cell
} // namespace walberla
