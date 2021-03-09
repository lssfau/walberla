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
//! \file Density.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename OutputType = float >
class DensityVTKWriter : public vtk::BlockCellDataWriter< OutputType >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   DensityVTKWriter( const ConstBlockDataID & pdf, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType >( id ), bdid_( pdf ), pdf_( nullptr ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( pdf_->getDensity(x,y,z) );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;

}; // class DensityVTKWriter



template< typename LatticeModel_T, typename OutputType = float >
class DensitySIVTKWriter : public vtk::BlockCellDataWriter< OutputType >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   DensitySIVTKWriter( const ConstBlockDataID & pdf, const real_t rho_SI, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType >( id ), bdid_( pdf ), pdf_( NULL ), rho_SI_( rho_SI ) {}

protected:

   void configure() { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( pdf_->getDensitySI( x, y, z, rho_SI_ ) );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;
   const real_t rho_SI_;

}; // class DensitySIVTKWriter



} // namespace lbm
} // namespace walberla
