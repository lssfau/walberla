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
//! \file DiffusionDirichlet.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
// @see Ginzbourg, I.   : Generic boundary conditions for lattice Boltzmann models and their application to advection and anisotropic disperion equations
// @see Like, L. et. al.: Boundary Conditions for thermal lattice Boltzmann equation method
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/ForceModel.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <vector>


namespace walberla {
namespace lbm {


template< typename LatticeModel_T, typename flag_t >
class DiffusionDirichlet : public Boundary<flag_t>
{
   static_assert( LatticeModel_T::compressible,                                                             "Only works with compressible models!" );
   //static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value), "Only works without additional forces!" );

   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   using ScalarField = GhostLayerField<real_t, 1>;

public:

   static const bool threadsafe = true;

   class ScalarConfiguration : public BoundaryConfiguration {
   public:
      virtual void val( real_t& _val, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const = 0;
   };

   class SingleScalarConfiguration : public ScalarConfiguration {
   public:
             SingleScalarConfiguration( const real_t _val = real_t(0) ) : val_( _val ) {}
      inline SingleScalarConfiguration( const Config::BlockHandle & config ){ val_ = ( config && config.isDefined( "val" ) ) ? config.getParameter< real_t >( "val" ) : real_c( 0. ); }

      const real_t & val() const { return val_; }

      void val( real_t& _val, cell_idx_t, cell_idx_t, cell_idx_t ) const override { _val = val(); }

      real_t & val(){ return val_; }

   private:
      real_t val_;
   };

   static shared_ptr<SingleScalarConfiguration> createConfiguration( const Config::BlockHandle & config ) { return make_shared<SingleScalarConfiguration>( config ); }

   inline DiffusionDirichlet( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, FlagField<flag_t> * const flagField = NULL );

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   inline void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & bc );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & bc );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & bc );

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

   inline real_t getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return sclField_->get(x,y,z); }

private:

   const FlagUID uid_;

   PDFField* const         pdfField_;
   shared_ptr<ScalarField> sclField_;

}; // class DiffusionDirichlet


template< typename LatticeModel_T, typename flag_t >
inline DiffusionDirichlet< LatticeModel_T, flag_t >::DiffusionDirichlet( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField, FlagField<flag_t> * const flagField ) :
   Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   if (flagField != nullptr)
      sclField_ = make_shared<ScalarField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), flagField->nrOfGhostLayers(), field::zyxf );
   else
      sclField_ = make_shared<ScalarField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField->nrOfGhostLayers(), field::zyxf );
}


template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << sclField_->get( x, y, z );
}


template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   buffer >> sclField_->get( x, y, z );
}


template< typename LatticeModel_T, typename flag_t >
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                const BoundaryConfiguration & bc )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const ScalarConfiguration * >( &bc ), &bc );
   WALBERLA_ASSERT_NOT_NULLPTR( sclField_ );

   const ScalarConfiguration & sclConfig = dynamic_cast< const ScalarConfiguration & >( bc );

   sclConfig.val( sclField_->get( x, y, z ), x, y, z );
}


template< typename LatticeModel_T, typename flag_t >
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::registerCells( const flag_t flag, const CellInterval & cells, const BoundaryConfiguration & bc )
{
   registerCells( flag, cells.begin(), cells.end(), bc );
}


template< typename LatticeModel_T, typename flag_t >
template< typename CellIterator >
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::registerCells( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                                 const BoundaryConfiguration & bc )
{
   for( auto cell = begin; cell != end; ++cell )
      registerCell( flag, cell->x(), cell->y(), cell->z(), bc );
}


template< typename LatticeModel_T, typename flag_t >
#ifndef NDEBUG
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                  const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void DiffusionDirichlet< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                  const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL  ( mask & this->mask_, this->mask_ ); 
   // only true if "this->mask_" only contains one single flag, which is the case for the current implementation of this boundary condition (DiffusionDirichlet)

   const real_t & bndVal = sclField_->get( nx, ny, nz );
   real_t &       target = pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) );

   target = - pdfField_->get( x, y, z, Stencil::idx[dir] ) + real_t(2) * bndVal * LatticeModel_T::w[Stencil::idx[dir]];
}


} // namespace lbm
} // namespace walberla
