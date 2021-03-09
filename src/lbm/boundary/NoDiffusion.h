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
//! \file NoDiffusion.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "NoSlip.h"


namespace walberla {
namespace lbm {

   // C++11 alias declaration not yet supported on all compilers
   //template< typename flag_t, typename Stencil > using NoDiffusion = NoSlip< flag_t, Stencil >;

   template< typename LatticeModel_T, typename flag_t > class NoDiffusion : public NoSlip< LatticeModel_T, flag_t >
   {
   protected:
      using PDFField = typename NoSlip<LatticeModel_T, flag_t>::PDFField;

   public:
      NoDiffusion( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField )
         : NoSlip< LatticeModel_T, flag_t >( boundaryUID, uid, pdfField ){}
   };

} // namespace lbm
} // namespace walberla
