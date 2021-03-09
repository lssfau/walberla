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
//! \file FreeDiffusion.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "FreeSlip.h"


namespace walberla {
namespace lbm {

   // C++11 alias declaration not yet supported on all compilers
   //template< typename FlagField_T, typename Stencil > using FreeDiffusion = FreeSlip< FlagField_T, Stencil >;

   template< typename LatticeModel_T, typename FlagField_T > class FreeDiffusion : public FreeSlip< LatticeModel_T, FlagField_T >
   {
   protected:
      using flag_t = typename FreeSlip<LatticeModel_T, FlagField_T>::flag_t;
      using PDFField = typename FreeSlip<LatticeModel_T, FlagField_T>::PDFField;

   public:
      FreeDiffusion( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField, const FlagField_T* const flagField, const flag_t domain )
         : FreeSlip< LatticeModel_T, FlagField_T >( boundaryUID, uid, pdfField, flagField, domain ){}
   };

} // namespace lbm
} // namespace walberla
