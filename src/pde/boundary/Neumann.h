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
//! \file Neumann.h
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
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

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"
#include "stencil/D3Q6.h"

#include <vector>
#include <limits>
#include <array>


namespace walberla {
namespace pde {



//**********************************************************************************************************************
/*!
*   \brief Class for setting Neumann domain boundaries for 1st and 2nd derivatives
*
*   ATTENTION: only works for cell-based setups where domain boundaries are located exactly in-between two cells !
*/
//**********************************************************************************************************************

template< typename PdeField >
class NeumannDomainBoundary
{
public:

   NeumannDomainBoundary( StructuredBlockStorage & blocks, const BlockDataID & fieldId ) :
      blocks_( blocks ), fieldId_( fieldId )
   {
      for( uint_t i = 0; i != stencil::D3Q6::Size; ++i )
      {
         includeBoundary_[i] = true;
         order_[i] = uint_t(1);
         value_[i] = real_t(0);
      }
      dx_[ stencil::D3Q6::idx[ stencil::W ] ] = blocks.dx();
      dx_[ stencil::D3Q6::idx[ stencil::E ] ] = blocks.dx();
      dx_[ stencil::D3Q6::idx[ stencil::S ] ] = blocks.dy();
      dx_[ stencil::D3Q6::idx[ stencil::N ] ] = blocks.dy();
      dx_[ stencil::D3Q6::idx[ stencil::B ] ] = blocks.dz();
      dx_[ stencil::D3Q6::idx[ stencil::T ] ] = blocks.dz();
   }

   void includeBoundary( const stencil::Direction & direction ) { includeBoundary_[stencil::D3Q6::idx[direction]] = true; }
   void excludeBoundary( const stencil::Direction & direction ) { includeBoundary_[stencil::D3Q6::idx[direction]] = false; }

   void setOrder( const uint_t order ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) order_[i] = order; }
   void setOrder( const stencil::Direction & direction, const uint_t order ) { order_[stencil::D3Q6::idx[direction]] = order; }

   void setValue( const real_t value ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) value_[i] = value; }
   void setValue( const stencil::Direction & direction, const real_t value ) { value_[stencil::D3Q6::idx[direction]] = value; }

   void setDx( const real_t dx ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) dx_[i] = dx; }
   void setDx( const stencil::Direction & direction, const real_t dx ) { dx_[stencil::D3Q6::idx[direction]] = dx; }

   void operator()();

protected:

   void apply( PdeField * p, const CellInterval & interval, const cell_idx_t cx, const cell_idx_t cy, const cell_idx_t cz,
               const uint_t order, const real_t value, const real_t dx ) const;



   StructuredBlockStorage & blocks_;
   BlockDataID fieldId_;

   bool includeBoundary_[ stencil::D3Q6::Size ];
   uint_t order_[ stencil::D3Q6::Size ];

   real_t value_[ stencil::D3Q6::Size ];
   real_t dx_[ stencil::D3Q6::Size ];

}; // class NeumannDomainBoundary



template< typename PdeField >
void NeumannDomainBoundary< PdeField >::operator()()
{
   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
   {
      PdeField * p = block->template getData< PdeField >( fieldId_ );

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::W ] ] && blocks_.atDomainXMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(-1), cell_idx_t(0), cell_idx_t(0),
                                 cell_idx_t(-1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(1), cell_idx_t(0), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::W ] ], value_[ stencil::D3Q6::idx[ stencil::W ] ], dx_[ stencil::D3Q6::idx[ stencil::W ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::E ] ] && blocks_.atDomainXMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_c(p->xSize()), cell_idx_t(0), cell_idx_t(0),
                                 cell_idx_c(p->xSize()), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(-1), cell_idx_t(0), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::E ] ], value_[ stencil::D3Q6::idx[ stencil::E ] ], dx_[ stencil::D3Q6::idx[ stencil::E ] ] );
      }

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::S ] ] && blocks_.atDomainYMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(-1), cell_idx_t(0),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_t(-1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(0), cell_idx_t(1), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::S ] ], value_[ stencil::D3Q6::idx[ stencil::S ] ], dx_[ stencil::D3Q6::idx[ stencil::S ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::N ] ] && blocks_.atDomainYMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_c(p->ySize()), cell_idx_t(0),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(0), cell_idx_t(-1), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::N ] ], value_[ stencil::D3Q6::idx[ stencil::N ] ], dx_[ stencil::D3Q6::idx[ stencil::N ] ] );
      }

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::B ] ] && blocks_.atDomainZMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(0), cell_idx_t(-1),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_t(-1) ),
                cell_idx_t(0), cell_idx_t(0), cell_idx_t(1),
                order_[ stencil::D3Q6::idx[ stencil::B ] ], value_[ stencil::D3Q6::idx[ stencil::B ] ], dx_[ stencil::D3Q6::idx[ stencil::B ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::T ] ] && blocks_.atDomainZMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(0), cell_idx_c(p->zSize()),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) ),
                cell_idx_t(0), cell_idx_t(0), cell_idx_t(-1),
                order_[ stencil::D3Q6::idx[ stencil::T ] ], value_[ stencil::D3Q6::idx[ stencil::T ] ], dx_[ stencil::D3Q6::idx[ stencil::T ] ] );
      }
   }
}



template< typename PdeField >
void NeumannDomainBoundary< PdeField >::apply( PdeField * p, const CellInterval & interval,
                                               const cell_idx_t cx, const cell_idx_t cy, const cell_idx_t cz,
                                               const uint_t order, const real_t value, const real_t dx ) const
{
   if( order == uint_t(1) )
   {
      if( isIdentical( value, real_t(0) ) )
      {
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,
            p->get(x,y,z) = p->get( x + cx, y + cy, z + cz );  // (dp / dx) == 0 _on_ the boundary
         )
      }
      else
      {
         const real_t vdx = value * dx;
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,
            p->get(x,y,z) = vdx + p->get( x + cx, y + cy, z + cz );  // (dp / dx) == value _on_ the boundary
         )
      }
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( order, uint_t(2) );

      if( isIdentical( value, real_t(0) ) )
      {      
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,

            const real_t pBoundary = p->get( x + cx, y + cy, z + cz );
            const real_t pInner    = p->get( x + cell_idx_t(2) * cx, y + cell_idx_t(2) * cy, z + cell_idx_t(2) * cz );

            const real_t boundaryValue = pBoundary + real_c(0.5) * ( pBoundary - pInner ); // extrapolation of value _on_ the boundary

            p->get(x,y,z) = real_t(2) * boundaryValue - pBoundary;  // (d^2 p / dx^2) == 0 _on_ the boundary
         )
      }
      else
      {
         const real_t vdx = value * real_c(0.25) * dx;
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,

            const real_t pBoundary = p->get( x + cx, y + cy, z + cz );
            const real_t pInner    = p->get( x + cell_idx_t(2) * cx, y + cell_idx_t(2) * cy, z + cell_idx_t(2) * cz );

            const real_t boundaryValue = pBoundary + real_c(0.5) * ( pBoundary - pInner ); // extrapolation of value _on_ the boundary

            p->get(x,y,z) = vdx + real_t(2) * boundaryValue - pBoundary;  // (d^2 p / dx^2) == value _on_ the boundary
         )
      }
   }
}



//**********************************************************************************************************************
/*!
*   \brief Neumann boundary handling for PDEs
*
*   This boundary condition imposes a Neumann condition with arbitrary values on a PDE.
*   It does so by modifying the right-hand side and the stencil field.
*   Anything that has internal copies of the stencil field (e.g. the multigrid V-cycle's coarsened stencils) is
*   responsible for updating its copies when boundary conditions are changed.
*
*   \tparam Stencil_T The stencil used for the discrete operator
*   \tparam flag_t The integer type used by the flag field
*/
//**********************************************************************************************************************
template< typename Stencil_T, typename flag_t >
class Neumann : public Boundary< flag_t >
{

   typedef GhostLayerField< real_t, 1 > Field_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

public:

   static const bool threadsafe = false;

   class NeumannBC : public BoundaryConfiguration {
   public:
             NeumannBC( const real_t & _neumannBC ) : neumannBC_( _neumannBC ) {}
      inline NeumannBC( const Config::BlockHandle & config );

      const real_t & neumannBC() const { return neumannBC_; }
      real_t & neumannBC() { return neumannBC_; }

   private:

      real_t neumannBC_;
   };

   static shared_ptr<NeumannBC> createConfiguration( const Config::BlockHandle & config ) { return make_shared<NeumannBC>( config ); }


  //*******************************************************************************************************************
  /*! Creates a Neumann boundary
   * \param boundaryUID the UID of the boundary condition
   * \param uid the UID of the flag that marks cells with this boundary condition
   * \param rhsField pointer to the right-hand side field, which will be adapted by this boundary condition
   * \param stencilField pointer to the operator stencil field. It should contain the stencil weights that don't take
   *                     into account the boundary conditions.
   * \param adaptBCStencilField pointer to the operator stencil field that will be adapted by this boundary condition. 
   *                            Initially, this field needs to contain the same values as \p stencilField.
   *                            This is the stencil field that should be passed to the actual PDE solver.
   * \param flagField pointer to the flag field
   * \param blocks
   *******************************************************************************************************************/
   inline Neumann( const BoundaryUID & boundaryUID, const FlagUID & uid, Field_T* const rhsField, const StencilField_T* const stencilField,
                     StencilField_T* const adaptBCStencilField, FlagField<flag_t> * const flagField, const StructuredBlockStorage& blocks );

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void afterBoundaryTreatment();

   template< typename Buffer_T >
   inline void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & neumannBC );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & neumannBC );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & neumannBC );

   void unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

   inline const real_t & getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const { return neumannBC_->get(x,y,z); }

private:

   const FlagUID uid_;
   const flag_t formerNeumann_;
   uint_t numDirtyCells_;

   std::array<real_t, Stencil_T::Q> dx_;
   Field_T* const                rhsField_;
   shared_ptr< Field_T >         neumannBC_;
   const StencilField_T * const  stencilField_;
   StencilField_T * const        adaptBCStencilField_;
   FlagField<flag_t> * const     flagField_;

}; // class Neumann



template< typename Stencil_T, typename flag_t >
inline Neumann< Stencil_T, flag_t >::NeumannBC::NeumannBC( const Config::BlockHandle & config  )
{
   neumannBC_ = ( config && config.isDefined( "val" ) ) ? config.getParameter<real_t>( "val" ) : real_c(0.0);
}



template< typename Stencil_T, typename flag_t >
inline Neumann< Stencil_T, flag_t >::Neumann( const BoundaryUID & boundaryUID, const FlagUID & uid, Field_T* const rhsField, const StencilField_T* const stencilField,
                                              StencilField_T* const adaptBCStencilField, FlagField<flag_t> * const flagField, const StructuredBlockStorage& blocks  )
                                              : Boundary<flag_t>( boundaryUID ), uid_( uid ), formerNeumann_ (flagField->getOrRegisterFlag("FormerNeumann")), numDirtyCells_(std::numeric_limits<uint_t>::max()),
                                                rhsField_( rhsField ), stencilField_ ( stencilField ), adaptBCStencilField_ ( adaptBCStencilField ), flagField_ (flagField)
{
   WALBERLA_ASSERT_NOT_NULLPTR( rhsField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( stencilField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( adaptBCStencilField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );

   WALBERLA_ASSERT_EQUAL( rhsField_->xyzSize(), stencilField_->xyzSize() );
#ifndef NDEBUG
   WALBERLA_FOR_ALL_CELLS_XYZ( stencilField_,
      for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
      {
         WALBERLA_ASSERT_IDENTICAL(stencilField_->get(x,y,z, dir.toIdx()), adaptBCStencilField_->get(x,y,z, dir.toIdx()));
      }
   )
#endif

   neumannBC_ = make_shared< Field_T >( rhsField_->xSize(), rhsField_->ySize(), rhsField_->zSize(), uint_t(1), field::fzyx );

   for(auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d ){
      dx_[d.toIdx()] = Vector3<real_t>(real_c(stencil::cx[d.toIdx()])*blocks.dx(), real_c(stencil::cy[d.toIdx()])*blocks.dy(), real_c(stencil::cz[d.toIdx()])*blocks.dz() ).sqrLength();
      // WALBERLA_LOG_DEVEL("dx in direction " << d.dirString() << ":" << dx_[d.toIdx()]);
   }

}


template< typename Stencil_T, typename flag_t >
void Neumann< Stencil_T, flag_t >::afterBoundaryTreatment() {

   if (numDirtyCells_>0 && numDirtyCells_ != std::numeric_limits<uint_t>::max()) {
      WALBERLA_LOG_WARNING("De-registering cells requires re-running Galerkin coarsening");
   }

   numDirtyCells_=0;
}

template< typename Stencil_T, typename flag_t >
template< typename Buffer_T >
inline void Neumann< Stencil_T, flag_t >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   buffer << neumannBC_->get( x, y, z );
}



template< typename Stencil_T, typename flag_t >
template< typename Buffer_T >
inline void Neumann< Stencil_T, flag_t >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   buffer >> neumannBC_->get( x, y, z );
   if( flagField_->isFlagSet( x, y, z, formerNeumann_ ) )
   {
      flagField_->removeFlag( x, y, z, formerNeumann_ );
      --numDirtyCells_;
   }
}



template< typename Stencil_T, typename flag_t >
inline void Neumann< Stencil_T, flag_t >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                       const BoundaryConfiguration & neumannBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const NeumannBC * >( &neumannBC ), &neumannBC );
   WALBERLA_ASSERT_NOT_NULLPTR( neumannBC_ );

   const NeumannBC & val = dynamic_cast< const NeumannBC & >( neumannBC );

   neumannBC_->get( x, y, z ) = val.neumannBC();

   if( flagField_->isFlagSet( x, y, z, formerNeumann_ ) )
   {
      flagField_->removeFlag( x, y, z, formerNeumann_ );
      --numDirtyCells_;
   }

}



template< typename Stencil_T, typename flag_t >
inline void Neumann< Stencil_T, flag_t >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & neumannBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const NeumannBC * >( &neumannBC ), &neumannBC );
   WALBERLA_ASSERT_NOT_NULLPTR( neumannBC_ );

   const NeumannBC & val = dynamic_cast< const NeumannBC & >( neumannBC );

   for( auto cell = neumannBC_->beginSliceXYZ( cells ); cell != neumannBC_->end(); ++cell ) {
      *cell = val.neumannBC();

      if( flagField_->isFlagSet( cell.cell(), formerNeumann_ ) )
      {
         flagField_->removeFlag( cell.cell(), formerNeumann_ );
         --numDirtyCells_;
      }
   }
}



template< typename Stencil_T, typename flag_t >
template< typename CellIterator >
inline void Neumann< Stencil_T, flag_t >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                                        const BoundaryConfiguration & neumannBC )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const NeumannBC * >( &neumannBC ), &neumannBC );
   WALBERLA_ASSERT_NOT_NULLPTR( neumannBC_ );

   const NeumannBC & val = dynamic_cast< const NeumannBC & >( neumannBC );

   for( auto cell = begin; cell != end; ++cell )
   {
      neumannBC_->get( cell->x(), cell->y(), cell->z() ) = val.neumannBC();

      if( flagField_->isFlagSet( cell->cell(), formerNeumann_ ) )
      {
         flagField_->removeFlag( cell->cell(), formerNeumann_ );
         --numDirtyCells_;
      }
   }

}

// Remark: This unregister function works only properly for D3Q7 stencils and convex domains!
template< typename Stencil_T, typename flag_t >
void Neumann< Stencil_T, flag_t >::unregisterCell( const flag_t, const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz )
{
   flagField_->addFlag( nx, ny, nz, formerNeumann_ );
   ++numDirtyCells_;

   Cell boundaryCell( nx, ny, nz );

   // Set stencil previously adapted to Neumann BC back to un-adapted state
   for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
   {
      Cell domainCell = boundaryCell - *d;
      if( adaptBCStencilField_->isInInnerPart( domainCell ) )
      {
         // restore original non-center stencil entry of neighboring non-boundary cell
         adaptBCStencilField_->get( domainCell, d.toIdx() ) = stencilField_->get( domainCell, d.toIdx() );

         // restore original center stencil entry of neighboring non-boundary cell
         adaptBCStencilField_->get( domainCell, Stencil_T::idx[stencil::C] ) -=  adaptBCStencilField_->get( domainCell, d.toIdx() );

      }
   }

}


//*************************************************************************************************
/*! \brief Treat one direction by adapting RHS according to Neumann boundary value.
 *
 */

template< typename Stencil_T, typename flag_t >
#ifndef NDEBUG
inline void Neumann< Stencil_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void Neumann< Stencil_T, flag_t >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{

   WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                             // current implementation of this boundary condition (Neumann)

   // WALBERLA_LOG_DEVEL("Neumann treatment at: " << Cell(x,y,z));

   // Adapt RHS to Neumann BC //
   rhsField_->get( x, y, z ) -= stencilField_->get( x, y, z, Stencil_T::idx[dir] ) * dx_[ Stencil_T::idx[dir] ] * neumannBC_->get( nx, ny, nz ); // possibly utilize that off-diagonal entries -1 anyway

   // Adapt Stencils to BCs (former adaptStencilsBC) //
   // Only required if any new BC cell was added or the BC type of any former BC cell has been changed
   if (numDirtyCells_>0) {

      // here not thread-safe! (workaround: mutex when adapting central entry)
      adaptBCStencilField_->get( x, y, z, Stencil_T::idx[stencil::C] ) += adaptBCStencilField_->get( x, y, z, Stencil_T::idx[dir] );
      adaptBCStencilField_->get( x, y, z, Stencil_T::idx[dir] ) = 0;

   }

}


} // namespace pde
} // namespace walberla
