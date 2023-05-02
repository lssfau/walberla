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
//! \file VelocityBoundary.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/FlagUID.h"
#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <vector>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename flag_t >
class VelocityBoundary : public Boundary<flag_t>
{
   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   using VelField = GhostLayerField<shared_ptr<std::vector<real_t>>, 1>;

public:

   static const bool threadsafe = true;

   class Configuration : public BoundaryConfiguration {
   public:
             Configuration( const Vector3< real_t > & _velocity, const real_t _density = real_t(1) ) :
                velocity_( _velocity ), density_( _density ) {}
             Configuration( const real_t _x, const real_t _y, const real_t _z, const real_t _density = real_t(1) ) :
                velocity_(_x,_y,_z), density_( _density ) {}
      inline Configuration( const Config::BlockHandle & config );

      const Vector3< real_t > & velocity() const { return velocity_; }
            real_t              density() const { return density_; }

      const real_t & x() const { return velocity_[0]; }
      const real_t & y() const { return velocity_[1]; }
      const real_t & z() const { return velocity_[2]; }

      Vector3< real_t > & velocity() { return velocity_; }
      real_t &            density() { return density_; }

      real_t & x() { return velocity_[0]; }
      real_t & y() { return velocity_[1]; }
      real_t & z() { return velocity_[2]; }

   private:

      Vector3< real_t > velocity_;
      real_t            density_;
   };

   static shared_ptr<Configuration> createConfiguration( const Config::BlockHandle & config ) { return make_shared<Configuration>( config ); }



   inline VelocityBoundary( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, FlagField<flag_t> * const flagField = NULL );

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

   inline const std::vector<real_t> & getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return *(vel_->get(x,y,z)); }

private:

   const FlagUID uid_;

   PDFField* const      pdfField_;
   shared_ptr<VelField> vel_;

}; // class VelocityBoundary



template< typename LatticeModel_T, typename flag_t >
inline VelocityBoundary< LatticeModel_T, flag_t >::Configuration::Configuration( const Config::BlockHandle & config  )
{
   velocity_[0] = ( config && config.isDefined( "x" ) ) ? config.getParameter<real_t>( "x" ) : real_c(0.0);
   velocity_[1] = ( config && config.isDefined( "y" ) ) ? config.getParameter<real_t>( "y" ) : real_c(0.0);
   velocity_[2] = ( config && config.isDefined( "z" ) ) ? config.getParameter<real_t>( "z" ) : real_c(0.0);
   density_ = ( config && config.isDefined( "rho" ) ) ? config.getParameter<real_t>( "rho" ) : real_t(1);
}



template< typename LatticeModel_T, typename flag_t >
inline VelocityBoundary< LatticeModel_T, flag_t >::VelocityBoundary( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField* const pdfField, FlagField<flag_t> * const flagField ) :

   Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   if (flagField != NULL)
      vel_ = make_shared<VelField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), flagField->nrOfGhostLayers(), field::fzyx );
   else
      vel_ = make_shared<VelField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField_->nrOfGhostLayers(), field::fzyx );
}



template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void VelocityBoundary< LatticeModel_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << *(vel_->get( x, y, z ));
}



template< typename LatticeModel_T, typename flag_t >
template< typename Buffer_T >
inline void VelocityBoundary< LatticeModel_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   std::vector< real_t > pdfs;
   buffer >> pdfs;
   vel_->get( x, y, z ) = walberla::make_shared< std::vector<real_t> >( pdfs );
}



template< typename LatticeModel_T, typename flag_t >
inline void VelocityBoundary< LatticeModel_T, flag_t >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                      const BoundaryConfiguration & velocity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &velocity ), &velocity );
   WALBERLA_ASSERT_NOT_NULLPTR( vel_ );

   const Configuration & vel = dynamic_cast< const Configuration & >( velocity );

   vel_->get( x, y, z ) = walberla::make_shared< std::vector<real_t> >( EquilibriumDistribution<LatticeModel_T>::get( vel.velocity(), vel.density() ) );
}



template< typename LatticeModel_T, typename flag_t >
inline void VelocityBoundary< LatticeModel_T, flag_t >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & velocity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &velocity ), &velocity );
   WALBERLA_ASSERT_NOT_NULLPTR( vel_ );

   const Configuration & vel = dynamic_cast< const Configuration & >( velocity );
   std::vector< real_t > pdfs = EquilibriumDistribution<LatticeModel_T>::get( vel.velocity(), vel.density() );

   for( auto cell = vel_->beginSliceXYZ( cells ); cell != vel_->end(); ++cell )
      *cell = walberla::make_shared< std::vector<real_t> >( pdfs );
}



template< typename LatticeModel_T, typename flag_t >
template< typename CellIterator >
inline void VelocityBoundary< LatticeModel_T, flag_t >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                       const BoundaryConfiguration & velocity )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &velocity ), &velocity );
   WALBERLA_ASSERT_NOT_NULLPTR( vel_ );

   const Configuration & vel = dynamic_cast< const Configuration & >( velocity );
   std::vector< real_t > pdfs = EquilibriumDistribution<LatticeModel_T>::get( vel.velocity(), vel.density() );

   for( auto cell = begin; cell != end; ++cell )
      vel_->get( cell->x(), cell->y(), cell->z() ) = walberla::make_shared< std::vector<real_t> >( pdfs );
}



template< typename LatticeModel_T, typename flag_t >
#ifndef NDEBUG
inline void VelocityBoundary< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                        const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void VelocityBoundary< LatticeModel_T, flag_t >::treatDirection( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const stencil::Direction dir,
                                                                        const cell_idx_t    nx, const cell_idx_t    ny, const cell_idx_t    nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                             // current implementation of this boundary condition (VelocityBoundary)

   WALBERLA_ASSERT( vel_->get( nx, ny, nz ).get() != nullptr );
   WALBERLA_ASSERT_LESS( Stencil::invDirIdx(dir), vel_->get( nx, ny, nz )->size() );

   pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = (*(vel_->get( nx, ny, nz )))[ Stencil::invDirIdx(dir) ];
}



} // namespace lbm
} // namespace walberla
