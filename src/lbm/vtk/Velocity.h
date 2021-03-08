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
//! \file Velocity.h
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
class VelocityVTKWriter : public vtk::BlockCellDataWriter< OutputType, 3 >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   VelocityVTKWriter( const ConstBlockDataID & pdfFieldId, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, 3 >( id ), bdid_( pdfFieldId ), pdf_( nullptr ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( (pdf_->getVelocity(x,y,z))[ uint_c(f) ] );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;

}; // class VelocityVTKWriter



template< typename LatticeModel_T, typename OutputType = float >
class VelocityMagnitudeVTKWriter : public vtk::BlockCellDataWriter< OutputType, 1 >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   VelocityMagnitudeVTKWriter( const ConstBlockDataID & pdfFieldId, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, 1 >( id ), bdid_( pdfFieldId ), pdf_( nullptr ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( pdf_->getVelocity(x,y,z).length() );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;

}; // class VelocityMagnitudeVTKWriter



template< typename LatticeModel_T, typename OutputType = float >
class VelocitySIVTKWriter : public vtk::BlockCellDataWriter< OutputType, 3 >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   VelocitySIVTKWriter( const ConstBlockDataID & pdfFieldId, const real_t dx_SI, const real_t dt_SI, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, 3 >( id ), bdid_( pdfFieldId ), pdf_( NULL ), dxDividedByDt_SI_( dx_SI / dt_SI ) {}

protected:

   void configure() { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( (pdf_->getVelocitySI( x, y, z, dxDividedByDt_SI_ ))[ uint_c(f) ] );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;
   const real_t dxDividedByDt_SI_;

}; // class VelocitySIVTKWriter



template< typename LatticeModel_T, typename OutputType = float >
class VelocitySIMagnitudeVTKWriter : public vtk::BlockCellDataWriter< OutputType, 1 >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   VelocitySIMagnitudeVTKWriter( const ConstBlockDataID & pdfFieldId, const real_t dx_SI, const real_t dt_SI, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, 1 >( id ), bdid_( pdfFieldId ), pdf_( NULL ), dxDividedByDt_SI_( dx_SI / dt_SI ) {}

protected:

   void configure() { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t /*f*/ )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );
      return numeric_cast< OutputType >( pdf_->getVelocitySI( x, y, z, dxDividedByDt_SI_ ).length() );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;
   const real_t dxDividedByDt_SI_;

}; // class VelocitySIMagnitudeVTKWriter



} // namespace lbm
} // namespace walberla
