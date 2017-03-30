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
//! \file TimeStep.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/NonUniformPackInfo.h"

namespace walberla {
namespace lbm {
namespace refinement {


class TimeStepPdfPackInfo : public blockforest::communication::NonUniformPackInfo
{
public:
   virtual bool optimizedEqualLevelCommunication() const = 0;
   virtual void optimizeEqualLevelCommunication( const bool value = true ) = 0;

   virtual bool optimizedForLinearExplosion() const = 0;
   virtual void optimizeForLinearExplosion( const bool value = true ) = 0;
};


} // namespace refinement
} // namespace lbm
} // namespace walberla
