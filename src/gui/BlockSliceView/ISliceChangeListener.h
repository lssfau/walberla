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
//! \file ISliceChangeListener.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/IBlock.h"


namespace walberla {
namespace gui {

   class ISliceChangeListener
   {
   public:
      virtual ~ISliceChangeListener() {}
      typedef int SliceID;
      virtual SliceID addSlice       ( IBlock * block, int sliceDim, double position ) = 0;
      virtual void    removeSlice    ( SliceID id )  = 0;
      virtual void    setSliceActive ( SliceID id, bool active = true ) = 0;
   };


} // namespace gui
} // namespace walberla


