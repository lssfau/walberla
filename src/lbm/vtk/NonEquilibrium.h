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
//! \file NonEquilibrium.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "core/debug/Debug.h"

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"

#include "vtk/BlockCellDataWriter.h"


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename OutputType = float >
class NonEqulibriumVTKWriter : public vtk::BlockCellDataWriter< OutputType, LatticeModel_T::Stencil::Size >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   NonEqulibriumVTKWriter( const ConstBlockDataID & pdf, const std::string & id ) :
      vtk::BlockCellDataWriter< OutputType, Stencil::Size >( id ), bdid_( pdf ), pdf_( nullptr ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );

      Vector3<real_t> v;

      real_t rho = pdf_->getDensityAndEquilibriumMomentumDensity( v, x, y, z );

      return numeric_cast< OutputType >( pdf_->get( x, y, z, f ) - lbm::EquilibriumDistribution< LatticeModel_T >::get( Stencil::dir[f], v, rho ) );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;

}; // class NonEqulibriumVTKWriter



} // namespace lbm
} // namespace walberla
