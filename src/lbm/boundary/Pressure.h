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
//! \file UBB.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
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
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <vector>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename flag_t >
class Pressure : public Boundary<flag_t>
{
   typedef PdfField< LatticeModel_T >        PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

   typedef GhostLayerField< real_t, 1 >  LatticeDensityField;

public:

   static const bool threadsafe = true;

   class LatticeDensity : public BoundaryConfiguration
   {
   public:
             LatticeDensity( const real_t & _latticeDensity ) : latticeDensity_( _latticeDensity ) {}
      inline LatticeDensity( const Config::BlockHandle & config );

      const real_t & latticeDensity() const { return latticeDensity_; }
            real_t   latticeDensity()       { return latticeDensity_; }
   private:
      real_t latticeDensity_;
   };

   static shared_ptr<LatticeDensity> createConfiguration( const Config::BlockHandle & config ) { return make_shared<LatticeDensity>( config ); }


   inline Pressure( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, FlagField<flag_t> * const flagField = NULL );

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   inline void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & velocity );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & velocity );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & velocity );

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

   inline real_t getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return latticeDensityField_->get(x,y,z); }

private:

   const FlagUID uid_;

   PDFField* const      pdfField_;
   shared_ptr<LatticeDensityField> latticeDensityField_;

}; // class Pressure



template< typename LatticeModel_T, typename flag_t >
inline Pressure< LatticeModel_T, flag_t >::LatticeDensity::LatticeDensity( const Config::BlockHandle & config  )
{
   latticeDensity_ =  ( config && config.isDefined( "latticeDensity" ) ) ? config.getParameter<real_t> ("latticeDensity") : real_t( 1.0 );
}



template< typename LatticeModel_T, typename flag_t >
inline Pressure< LatticeModel_T, flag_t>::Pressure( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, FlagField<flag_t> * const flagField ) :

   Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   if (flagField != NULL)
      latticeDensityField_ = make_shared<LatticeDensityField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), flagField->nrOfGhostLayers(), field::zyxf );
   else
      latticeDensityField_ = make_shared<LatticeDensityField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField_->nrOfGhostLayers(), field::zyxf );
}



template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void Pressure< LatticeModel_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << latticeDensityField_->get( x, y, z );
}



template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void Pressure< LatticeModel_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   buffer >> latticeDensityField_->get( x, y, z );
}



template< typename LatticeModel_T, typename flag_t >
inline void Pressure< LatticeModel_T, flag_t >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                            const BoundaryConfiguration & latticeDensity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const LatticeDensity * >( &latticeDensity ), &latticeDensity );
   WALBERLA_ASSERT_NOT_NULLPTR( latticeDensityField_ );

   const LatticeDensity & rhoConfig = dynamic_cast< const LatticeDensity & >( latticeDensity );

   latticeDensityField_->get( x, y, z ) = rhoConfig.latticeDensity();
}



template< typename LatticeModel_T, typename flag_t >
inline void Pressure< LatticeModel_T, flag_t >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & latticeDensity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const LatticeDensity * >( &latticeDensity ), &latticeDensity );
   WALBERLA_ASSERT_NOT_NULLPTR( latticeDensityField_ );

   const LatticeDensity & rhoConfig = dynamic_cast< const LatticeDensity & >( latticeDensity );

   for( auto cell = latticeDensityField_->beginSliceXYZ( cells ); cell != latticeDensityField_->end(); ++cell )
      *cell = rhoConfig.latticeDensity();
}



template< typename LatticeModel_T, typename flag_t >
template< typename CellIterator >
inline void Pressure< LatticeModel_T, flag_t >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                                             const BoundaryConfiguration & latticeDensity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const LatticeDensity * >( &latticeDensity ), &latticeDensity );
   WALBERLA_ASSERT_NOT_NULLPTR( latticeDensityField_ );

   const LatticeDensity & rhoConfig = dynamic_cast< const LatticeDensity & >( latticeDensity );

   for( auto cell = begin; cell != end; ++cell )
      latticeDensityField_->get( cell->x(), cell->y(), cell->z() ) = rhoConfig.latticeDensity();
}



template< typename LatticeModel_T, typename flag_t >
#ifndef NDEBUG
inline void Pressure< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void Pressure< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );

   WALBERLA_ASSERT_UNEQUAL( ( mask & this->mask_ ), numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL( ( mask & this->mask_ ), this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                 // current implementation of this boundary condition (Pressure)
   Vector3<real_t> u = pdfField_->getVelocity(x,y,z);

   // result will be streamed to (x,y,z, stencil::inverseDir[d]) during sweep
   pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) =
      - pdfField_->get( x, y, z, Stencil::idx[dir] )                   //anti-bounce-back
      + real_t(2) * EquilibriumDistribution<LatticeModel_T>::getSymmetricPart( dir, u,  latticeDensityField_->get(nx,ny,nz) ); //pressure term
}



} // namespace lbm
} // namespace walberla
