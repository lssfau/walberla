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
//! \file PdfFieldSyncPackInfo.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "field/refinement/PackInfo.h"


namespace walberla {
namespace lbm {
namespace refinement {

template <typename LatticeModel_T>
class PdfFieldSyncPackInfo : public field::refinement::PackInfo< PdfField<LatticeModel_T>, typename LatticeModel_T::Stencil >
{
public:
   PdfFieldSyncPackInfo( const BlockDataID & fieldId )
   : field::refinement::PackInfo< PdfField<LatticeModel_T>, typename LatticeModel_T::Stencil >( fieldId )
   {
   }
};

} // namespace refinement
} // namespace lbm
} // namespace walberla
