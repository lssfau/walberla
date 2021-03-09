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
//! \file PressureTensor.h
//! \ingroup lbm
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//================

#pragma once

#include "lbm/field/PdfField.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"
#include "lbm/field/PressureTensor.h"


namespace walberla {
namespace lbm {


template< typename LatticeModel_T, typename OutputType = float >
class PressureTensorVTKWriter : public vtk::BlockCellDataWriter< OutputType, 9 >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   PressureTensorVTKWriter( const ConstBlockDataID & pdfFieldId, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, 9 >( id ), bdid_( pdfFieldId ), pdf_( NULL ) {}

protected:

   void configure() { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( (pdf_->getPressureTensor(x,y,z))[ uint_c(f) ] );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;

}; // class PressureTensorVTKWriter

} // namespace lbm
} // namespace walberla
