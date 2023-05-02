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
//! \file Curved.h
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

#include <array>
#include <vector>

namespace walberla {
namespace lbm {



// delta = "distance to wall" / "distance between cell centers"

template< typename LatticeModel_T, typename FlagField_T >
class Curved : public Boundary< typename FlagField_T::flag_t >
{
   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagField_T::flag_t;

   using WeightField = GhostLayerField<shared_ptr<std::array<real_t, Stencil::Size>>, 1>;

public:

   static const bool threadsafe = true;

   class Configuration : public BoundaryConfiguration {
   public:
             Configuration( const std::vector< real_t > & _deltas ) : deltas_( _deltas ) { WALBERLA_ASSERT_EQUAL( _deltas.size(), Stencil::Size ); }
      inline Configuration( const Config::BlockHandle & config );

      const std::vector< real_t > & deltas() const { return deltas_; }
            std::vector< real_t > & deltas()       { return deltas_; }

   private:

      std::vector< real_t > deltas_;
   };

   static shared_ptr<Configuration> createConfiguration( const Config::BlockHandle & config ) { return make_shared<Configuration>( config ); }


   
   inline Curved( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                  const FlagField_T * const flagField, const flag_t domain );

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   inline void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & deltas );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & deltas );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & deltas );

   inline void unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

   inline const typename WeightField::value_type & getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return weights_->get(x,y,z); }

private:

   inline static real_t deltaToWeight( const real_t delta ) { return ( real_t( 1 ) - real_t( 2 ) * delta ) / ( real_t( 1 ) + real_t( 2 ) * delta ); }
   inline static shared_ptr< std::array<real_t, Stencil::Size> > makeDeltaArray()
   {
      return walberla::shared_ptr< std::array<real_t, Stencil::Size> >( new std::array<real_t, Stencil::Size> );
   }

   const FlagUID uid_;
   
   flag_t domainMask_;

         PDFField *    const pdfField_;
   const FlagField_T * const flagField_;   
      
   shared_ptr<WeightField> weights_; 

}; // class Curved



template< typename LatticeModel_T, typename FlagField_T >
inline Curved< LatticeModel_T, FlagField_T >::Configuration::Configuration( const Config::BlockHandle & config  )
{
   deltas_.resize( Stencil::Size, real_c(0.5) );
   for( uint_t i = 0; i < Stencil::Size; ++i )
   {
      std::ostringstream oss;
      oss << "delta" << i;
      deltas_[i] = ( config && config.isDefined( oss.str() ) ) ? config.getParameter<real_t>( oss.str() ) : real_c(0.5);
   }
}



template< typename LatticeModel_T, typename FlagField_T >
inline Curved< LatticeModel_T, FlagField_T >::Curved( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                                                      const FlagField_T * const flagField, const flag_t domain ) :

   Boundary<flag_t>( boundaryUID ), uid_( uid ), domainMask_(domain), pdfField_( pdfField ), flagField_( flagField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
   WALBERLA_ASSERT( flagField_->isRegistered( domainMask_ )  );   
   weights_ = make_shared<WeightField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), flagField_->nrOfGhostLayers(), field::fzyx );
}



template< typename LatticeModel_T, typename FlagField_T >
template< typename Buffer_T >
inline void Curved< LatticeModel_T, FlagField_T >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << *(weights_->get( x, y, z ));
}



template< typename LatticeModel_T, typename FlagField_T >
template< typename Buffer_T >
inline void Curved< LatticeModel_T, FlagField_T >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   weights_->get( x, y, z ) = makeDeltaArray();
   buffer >> *( weights_->get( x, y, z ) );
}



template< typename LatticeModel_T, typename FlagField_T >
inline void Curved< LatticeModel_T, FlagField_T >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                 const BoundaryConfiguration & deltas )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &deltas ), &deltas );
   WALBERLA_ASSERT_NOT_NULLPTR( weights_ );

   const Configuration & config = dynamic_cast< const Configuration & >( deltas );

   weights_->get( x, y, z ) = makeDeltaArray();
   
   std::transform( config.deltas().begin(), config.deltas().end(), weights_->get( x, y, z )->begin(), deltaToWeight );
}



template< typename LatticeModel_T, typename FlagField_T >
inline void Curved< LatticeModel_T, FlagField_T >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & deltas )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &deltas ), &deltas );
   WALBERLA_ASSERT_NOT_NULLPTR( weights_ );

   const Configuration & config = dynamic_cast< const Configuration & >( deltas );

   for( auto cell = weights_->beginSliceXYZ( cells ); cell != weights_->end(); ++cell )
   {
      *cell = makeDeltaArray();
      std::transform( config.deltas().begin(), config.deltas().end(), ( *cell )->begin(), deltaToWeight );
   }
}



template< typename LatticeModel_T, typename FlagField_T >
template< typename CellIterator >
inline void Curved< LatticeModel_T, FlagField_T >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                  const BoundaryConfiguration & deltas )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Configuration * >( &deltas ), &deltas );
   WALBERLA_ASSERT_NOT_NULLPTR( weights_ );

   const Configuration & config = dynamic_cast< const Configuration & >( deltas );

   for( auto cell = begin; cell != end; ++cell )
   {
      weights_->get( cell->x(), cell->y(), cell->z() ) = makeDeltaArray();
      std::transform( config.deltas().begin(), config.deltas().end(), weights_->get( cell->x(), cell->y(), cell->z() )->begin(), deltaToWeight );
   }
}



template< typename LatticeModel_T, typename FlagField_T >
inline void Curved< LatticeModel_T, FlagField_T >::unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   weights_->get( x, y, z ).reset();
}



template< typename LatticeModel_T, typename FlagField_T >
#ifndef NDEBUG
inline void Curved< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                   const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void Curved< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                   const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                             // current implementation of this boundary condition (Curved)
   
   const cell_idx_t xff = x - cell_idx_c( stencil::cx[ dir ] );
   const cell_idx_t yff = y - cell_idx_c( stencil::cy[ dir ] );
   const cell_idx_t zff = z - cell_idx_c( stencil::cz[ dir ] );
   
   if( !flagField_->isPartOfMaskSet( xff, yff, zff, domainMask_ ) )
   {
      // default no slip
      pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] );
   }
   else
   {
      WALBERLA_ASSERT( weights_->get( nx, ny, nz ).get() != nullptr );
      WALBERLA_ASSERT_LESS( Stencil::invDirIdx(dir), weights_->get( nx, ny, nz )->size() );
      
      // linear multi reflection model without non-equilibrium

      const real_t weight = ( *( weights_->get( nx, ny, nz ) ) )[Stencil::invDirIdx( dir )];

      pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] )
                                                              + weight * ( pdfField_->get( xff, yff, zff, Stencil::idx[dir] )
                                                                           - pdfField_->get( x, y, z, Stencil::invDirIdx(dir) ) );      
   }
}


} // namespace lbm
} // namespace walberla
