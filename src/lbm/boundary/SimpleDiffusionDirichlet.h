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
//! \file SimpleDiffusionDirichlet.h
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
class SimpleDiffusionDirichlet : public Boundary<flag_t>
{
   static_assert( LatticeModel_T::compressible,                                                             "Only works with compressible models!" );
   //static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value), "Only works without additional forces!" );

   typedef PdfField< LatticeModel_T >        PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

public:

   static const bool threadsafe = true;

   class ScalarConfiguration : public BoundaryConfiguration {
   public:
             ScalarConfiguration( const real_t _val = real_t(0) ) : val_( _val ) {}
      inline ScalarConfiguration( const Config::BlockHandle & config ){ val_ = ( config && config.isDefined( "val" ) ) ? config.getParameter< real_t >( "val" ) : real_t(0); }

      const real_t & val() const { return val_; }
            real_t & val()       { return val_; }

   private:
      real_t val_;
   };

   static shared_ptr<ScalarConfiguration> createConfiguration( const Config::BlockHandle & config ) { return make_shared<ScalarConfiguration>( config ); }

   inline SimpleDiffusionDirichlet( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, const real_t val = real_t(0) );

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

private:

   const FlagUID uid_;

   PDFField* const pdfField_;

   real_t val_;
   bool   init_;

}; // class SimpleDiffusionDirichlet

template< typename LatticeModel_T, typename flag_t >
inline SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::SimpleDiffusionDirichlet( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField, const real_t val ) :
   Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), val_( val ), init_(false) 
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ ); 
}


template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const
{
   buffer << val_;
}


template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t )
{
   real_t tmp;
   buffer >> tmp;

   if( init_ )
      WALBERLA_ASSERT_FLOAT_EQUAL( val_, tmp );

   init_ = true;
   val_  = tmp;
}


template< typename LatticeModel_T, typename flag_t >
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                                                                                const BoundaryConfiguration & bc )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const ScalarConfiguration * >( &bc ), &bc );

   const ScalarConfiguration & sclConfig = dynamic_cast< const ScalarConfiguration & >( bc );

   if( init_ )
      WALBERLA_ASSERT_FLOAT_EQUAL( val_, sclConfig.val() );

   init_ = true;
   val_  = sclConfig.val();
}


template< typename LatticeModel_T, typename flag_t >
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::registerCells( const flag_t flag, const CellInterval & cells, const BoundaryConfiguration & bc )
{
   registerCells( flag, cells.begin(), cells.end(), bc );
}


template< typename LatticeModel_T, typename flag_t >
template< typename CellIterator >
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::registerCells( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                                 const BoundaryConfiguration & bc )
{
   for( auto cell = begin; cell != end; ++cell )
      registerCell( flag, cell->x(), cell->y(), cell->z(), bc );
}


template< typename LatticeModel_T, typename flag_t >
#ifndef NDEBUG
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                  const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void SimpleDiffusionDirichlet< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                  const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL  ( mask & this->mask_, this->mask_ ); 
   // only true if "this->mask_" only contains one single flag, which is the case for the current implementation of this boundary condition (SimpleDiffusionDirichlet)

   pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = real_t(2) * val_ * LatticeModel_T::w[Stencil::idx[dir]] - pdfField_->get( x, y, z, Stencil::idx[dir] );
}


} // namespace lbm
} // namespace walberla
