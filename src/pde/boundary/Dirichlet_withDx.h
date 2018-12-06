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
//! \file Dirichlet.h
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "field/Field.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"

#include <vector>
#include <array>


namespace walberla {
namespace pde {



template< typename Stencil_T, typename flag_t >
class Dirichlet : public Boundary< flag_t >
{

   typedef GhostLayerField< real_t, 1 > Field_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

public:

   static const bool threadsafe = false;

   class DirichletBC : public BoundaryConfiguration {
   public:
             DirichletBC( const real_t & _dirichletBC ) : dirichletBC_( _dirichletBC ) {}
      inline DirichletBC( const Config::BlockHandle & config );

      const real_t & dirichletBC() const { return dirichletBC_; }
      real_t & dirichletBC() { return dirichletBC_; }

   private:

      real_t dirichletBC_;
   };

   static shared_ptr<DirichletBC> createConfiguration( const Config::BlockHandle & config ) { return make_shared<DirichletBC>( config ); }



   inline Dirichlet( const BoundaryUID & boundaryUID, const FlagUID & uid, Field_T* const rhsField, const StencilField_T* const stencilField,
                     StencilField_T* const adaptBCStencilField, FlagField<flag_t> * const flagField, const StructuredBlockStorage& blocks );

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void afterBoundaryTreatment();

   template< typename Buffer_T >
   inline void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & dirichletBC );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & dirichletBC );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & dirichletBC );

   void unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

   inline const real_t & getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return dirichletBC_->get(x,y,z); }

private:

   const FlagUID uid_;
   const flag_t formerDirichlet_;
   uint_t numDirtyCells_;

   std::array<real_t, Stencil_T::Q> dx_;
   Field_T* const                rhsField_;
   shared_ptr< Field_T >         dirichletBC_;
   const StencilField_T * const  stencilField_;
   StencilField_T * const        adaptBCStencilField_;
   FlagField<flag_t> * const     flagField_;

}; // class Dirichlet



template< typename Stencil_T, typename flag_t >
inline Dirichlet< Stencil_T, flag_t >::DirichletBC::DirichletBC( const Config::BlockHandle & config  )
{
   dirichletBC_ = ( config && config.isDefined( "val" ) ) ? config.getParameter<real_t>( "val" ) : real_c(0.0);
}



template< typename Stencil_T, typename flag_t >
inline Dirichlet< Stencil_T, flag_t >::Dirichlet( const BoundaryUID & boundaryUID, const FlagUID & uid, Field_T* const rhsField, const StencilField_T* const stencilField, StencilField_T* const adaptBCStencilField, FlagField<flag_t> * const flagField, const StructuredBlockStorage& blocks ) :
   Boundary<flag_t>( boundaryUID ), uid_( uid ), formerDirichlet_ (flagField->getOrRegisterFlag("FormerDirichlet")), numDirtyCells_(0), rhsField_( rhsField ), stencilField_ ( stencilField ), adaptBCStencilField_ ( adaptBCStencilField ), flagField_ (flagField)
{
   WALBERLA_ASSERT_NOT_NULLPTR( rhsField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( stencilField_ );

   WALBERLA_ASSERT_EQUAL( rhsField_->xyzSize(), stencilField_->xyzSize() );

   dirichletBC_ = make_shared< Field_T >( rhsField_->xSize(), rhsField_->ySize(), rhsField_->zSize(), uint_t(1), field::zyxf );

   for(auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d ){
      dx_[d.toIdx()] = Vector3<real_t>(stencil::cx[d.toIdx()]*blocks.dx(), stencil::cy[d.toIdx()]*blocks.dy(), stencil::cz[d.toIdx()]*blocks.dz() ).sqrLength();
      WALBERLA_LOG_DEVEL("dx in direction " << d.dirString() << ":" << dx_[d.toIdx()]);
   }


}


template< typename Stencil_T, typename flag_t >
void Dirichlet< Stencil_T, flag_t >::afterBoundaryTreatment() {

   if (numDirtyCells_>0) {
      WALBERLA_LOG_WARNING("De-registering cells requires re-running Galerkin coarsening");
   }

   numDirtyCells_=0;
}

template< typename Stencil_T, typename flag_t >
template< typename Buffer_T >
inline void Dirichlet< Stencil_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << dirichletBC_->get( x, y, z );
}



template< typename Stencil_T, typename flag_t >
template< typename Buffer_T >
inline void Dirichlet< Stencil_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   buffer >> dirichletBC_->get( x, y, z );
   if( flagField_->isFlagSet( x, y, z, formerDirichlet_ ) )
   {
      flagField_->removeFlag( x, y, z, formerDirichlet_ );
      --numDirtyCells_;
   }
}



template< typename Stencil_T, typename flag_t >
inline void Dirichlet< Stencil_T, flag_t >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                       const BoundaryConfiguration & dirichletBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const DirichletBC * >( &dirichletBC ), &dirichletBC );
   WALBERLA_ASSERT_NOT_NULLPTR( dirichletBC_ );

   const DirichletBC & val = dynamic_cast< const DirichletBC & >( dirichletBC );

   dirichletBC_->get( x, y, z ) = val.dirichletBC();

   if( flagField_->isFlagSet( x, y, z, formerDirichlet_ ) )
   {
      flagField_->removeFlag( x, y, z, formerDirichlet_ );
      --numDirtyCells_;
   }

}



template< typename Stencil_T, typename flag_t >
inline void Dirichlet< Stencil_T, flag_t >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & dirichletBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const DirichletBC * >( &dirichletBC ), &dirichletBC );
   WALBERLA_ASSERT_NOT_NULLPTR( dirichletBC_ );

   const DirichletBC & val = dynamic_cast< const DirichletBC & >( dirichletBC );

   for( auto cell = dirichletBC_->beginSliceXYZ( cells ); cell != dirichletBC_->end(); ++cell ) {
      *cell = val.dirichletBC();

      if( flagField_->isFlagSet( cell.cell(), formerDirichlet_ ) )
      {
         flagField_->removeFlag( cell.cell(), formerDirichlet_ );
         --numDirtyCells_;
      }
   }
}



template< typename Stencil_T, typename flag_t >
template< typename CellIterator >
inline void Dirichlet< Stencil_T, flag_t >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                                        const BoundaryConfiguration & dirichletBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const DirichletBC * >( &dirichletBC ), &dirichletBC );
   WALBERLA_ASSERT_NOT_NULLPTR( dirichletBC_ );

   const DirichletBC & val = dynamic_cast< const DirichletBC & >( dirichletBC );

   for( auto cell = begin; cell != end; ++cell )
   {
      dirichletBC_->get( cell->x(), cell->y(), cell->z() ) = val.dirichletBC();

      if( flagField_->isFlagSet( cell->cell(), formerDirichlet_ ) )
      {
         flagField_->removeFlag( cell->cell(), formerDirichlet_ );
         --numDirtyCells_;
      }
   }

}


template< typename Stencil_T, typename flag_t >
void Dirichlet< Stencil_T, flag_t >::unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   flagField_->addFlag( x,y,z, formerDirichlet_ );
   ++numDirtyCells_;

   // Set stencil adapted to BCs back to unadapted state
   for(auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d ){
      adaptBCStencilField_->get(x,y,z,d.toIdx()) = stencilField_->get(x,y,z,d.toIdx());
   }

}


//*************************************************************************************************
/*! \brief Treat one direction by adapting RHS according to Dirichlet boundary value.
 *
 */

template< typename Stencil_T, typename flag_t >
#ifndef NDEBUG
inline void Dirichlet< Stencil_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void Dirichlet< Stencil_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{

   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                             // current implementation of this boundary condition (Dirichlet)

   // Adapt RHS to Dirichlet BC //
   rhsField_->get( x, y, z ) -= stencilField_->get( x, y, z, Stencil_T::idx[dir] ) * real_t(2) * dx_[ Stencil_T::idx[dir] ] * dirichletBC_->get( nx, ny, nz ); // possibly utilize that off-diagonal entries -1 anyway

   // WALBERLA_LOG_DEVEL("Adapt RHS to Dirichlet value " << dirichletBC_->get( nx, ny, nz ) << " on cell " << Cell(x,y,z) << " for stencil entry " << stencilField_->get( x, y, z, Stencil_T::idx[dir] ) );

   // Adapt Stencils to BCs (former adaptStencilsBC) //
   // Only required if any new BC cell was added or the BC type of any former BC cell has been changed
   if (numDirtyCells_>0) {

      // here not thread-safe! (workaround: mutex when adapting central entry)
      adaptBCStencilField_->get( x, y, z, Stencil_T::idx[stencil::C] ) -= adaptBCStencilField_->get( x, y, z, Stencil_T::idx[dir] );
      adaptBCStencilField_->get( x, y, z, Stencil_T::idx[dir] ) = 0;

   }

}



} // namespace pde
} // namespace walberla
