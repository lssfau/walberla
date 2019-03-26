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
//! \file NoSlip.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
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



template< typename LatticeModel_T, typename flag_t, bool StoreForce = false >
class NoSlip : public Boundary<flag_t>
{
protected:

   typedef PdfField< LatticeModel_T >        PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

   typedef GhostLayerField< Vector3<real_t>, 1 > ForceField;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   NoSlip( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField ) :
      Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );

         if (StoreForce)
            force_ = make_shared<ForceField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField_->nrOfGhostLayers(), field::zyxf );
      }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const
   {
      if (StoreForce)
         force_->setWithGhostLayer( Vector3<real_t>() );
   }
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
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (NoSlip)

      const real_t pdf_old = pdfField_->get( x, y, z, Stencil::idx[dir] );

      pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] );

      if (StoreForce && pdfField_->isInInnerPart( Cell(x,y,z) ))
      {
         const real_t forceMEM = pdf_old + pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) );
         Vector3<real_t> force( real_c( stencil::cx[dir] ) * forceMEM,
                                real_c( stencil::cy[dir] ) * forceMEM,
                                real_c( stencil::cz[dir] ) * forceMEM );
         force_->get( nx, ny, nz ) += force;
      }
   }

   const typename ForceField::value_type & getForce( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      static_assert(StoreForce, "this member function is only available if the third template argument on the class is true");
      return force_->get(x,y,z);
   }

private:

   const FlagUID uid_;

   PDFField* const pdfField_;

   shared_ptr<ForceField> force_;

}; // class NoSlip



} // namespace lbm
} // namespace walberla
