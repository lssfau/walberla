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
//! \file Outlet.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Ehsan Fattahi <ehsan.fattahi@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"

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



template< typename LatticeModel_T, typename FlagField_T, int Numerator = 2, int Denominator = 1 >
class Outlet : public Boundary< typename FlagField_T::flag_t > {

protected:
   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagField_T::flag_t;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   Outlet( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField, const FlagField_T* const flagField, const flag_t domain ) :
      Boundary<flag_t>( boundaryUID ), uid_( uid ), domainMask_(domain), pdfField_( pdfField ), flagField_( flagField )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
      WALBERLA_ASSERT( flagField_->isRegistered( domainMask_ )  );
   }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   void registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration& ) {}
   void registerCells( const flag_t, const CellInterval&, const BoundaryConfiguration& ) {}
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator&, const CellIterator&, const BoundaryConfiguration& ) {}

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
   inline void treatDirection( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const stencil::Direction dir,
                               const cell_idx_t nx,    const cell_idx_t ny,    const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (Outlet)

      // the following computes the normal of the outlet boundary (by searching for fluid neighbors)

      // inverse direction of 'dir' as lattice vector

      const int ix = stencil::cx[ stencil::inverseDir[dir] ];
      const int iy = stencil::cy[ stencil::inverseDir[dir] ];
      const int iz = stencil::cz[ stencil::inverseDir[dir] ];

      int wnx = 0; // compute "normal" vector of outlet wall
      int wny = 0;
      int wnz = 0;

      WALBERLA_ASSERT( flagField_->xyzSizeWithGhostLayer().contains( nx+ix, ny, nz ) );
      WALBERLA_ASSERT( flagField_->xyzSizeWithGhostLayer().contains( nx, ny+iy, nz ) );
      WALBERLA_ASSERT( flagField_->xyzSizeWithGhostLayer().contains( nx, ny, nz+iz ) );

      if( flagField_->isPartOfMaskSet( nx+ix, ny, nz, domainMask_ ) ) wnx = ix;
      if( flagField_->isPartOfMaskSet( nx, ny+iy, nz, domainMask_ ) ) wny = iy;
      if( flagField_->isPartOfMaskSet( nx, ny, nz+iz, domainMask_ ) ) wnz = iz;

      // concave corner (neighbors are non-fluid)
      if( wnx == 0 && wny == 0 && wnz == 0 )
      {
         wnx = ix;
         wny = iy;
         wnz = iz;
      }

      int wnx2 = wnx;
      int wny2 = wny;
      int wnz2 = wnz;

      WALBERLA_ASSERT( flagField_->xyzSizeWithGhostLayer().contains( nx+2*wnx, ny+2*wny, nz+2*wnz ) );

      if( flagField_->isPartOfMaskSet( nx+2*wnx, ny+2*wny, nz+2*wnz, domainMask_ ) )
      {
         wnx2 = 2 * wnx;
         wny2 = 2 * wny;
         wnz2 = 2 * wnz;
      }

      WALBERLA_ASSERT( pdfField_->xyzSizeWithGhostLayer().contains( nx + wnx,  ny + wny,  nz + wnz  ) );
      WALBERLA_ASSERT( pdfField_->xyzSizeWithGhostLayer().contains( nx + wnx2, ny + wny2, nz + wnz2 ) );
      
      if( Denominator == 1 )
      {
         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) =
            real_t( Numerator ) * pdfField_->get( nx + wnx,  ny + wny,  nz + wnz,  Stencil::invDirIdx(dir) )
                                - pdfField_->get( nx + wnx2, ny + wny2, nz + wnz2, Stencil::invDirIdx(dir) );
      }
      else
      {
         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) =
            ( real_t( Numerator ) * pdfField_->get( nx + wnx,  ny + wny,  nz + wnz,  Stencil::invDirIdx(dir) )
                                  - pdfField_->get( nx + wnx2, ny + wny2, nz + wnz2, Stencil::invDirIdx(dir) ) ) / real_t( Denominator );
      }
   }

private:

   const FlagUID uid_;

   flag_t domainMask_;

         PDFField *    const pdfField_;
   const FlagField_T * const flagField_;

}; // class Outlet



} // namespace lbm
} // namespace walberla
