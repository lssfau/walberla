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
//! \file SimplePressure.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/ForceModel.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <vector>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename flag_t >
class SimplePressure : public Boundary<flag_t>
{
   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   using Vec3Real = Vector3<real_t>;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   SimplePressure( const BoundaryUID& boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                   const real_t latticeDensity )
      : Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfs_( pdfField ), latticeDensity_( latticeDensity )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfs_  );
   }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   void registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration& ) {}
   void registerCells( const flag_t, const CellInterval&, const BoundaryConfiguration& ) const {}
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator&, const CellIterator&, const BoundaryConfiguration& ) const {}

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   void setLatticeDensity( const real_t newLatticeDensity ) {
      latticeDensity_ = newLatticeDensity;
   }

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );

      WALBERLA_ASSERT_UNEQUAL( ( mask & this->mask_ ), numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( ( mask & this->mask_ ), this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                    // current implementation of this boundary condition (SimplePressure)
      Vec3Real u = pdfs_->getVelocity(x,y,z);

      // result will be streamed to (x,y,z, stencil::inverseDir[d]) during sweep
      pdfs_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) =
         - pdfs_->get( x, y, z, Stencil::idx[dir] )                   //anti-bounce-back
         + real_t(2) * EquilibriumDistribution<LatticeModel_T>::getSymmetricPart( dir, u, latticeDensity_ ); //pressure term
   }

protected:

   FlagUID uid_;

   PDFField * pdfs_;
   real_t latticeDensity_;

}; // class SimplePressure



} // namespace lbm
} // namespace walberla
