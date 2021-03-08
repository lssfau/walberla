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
//! \file SimplePAB.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/ForceModel.h"

#include "boundary/Boundary.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <type_traits>

#include <vector>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename FlagFieldT >
class SimplePAB : public Boundary<typename FlagFieldT::flag_t>
{
   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );

   using PDFField = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagFieldT::flag_t;

   using Vec3Real = Vector3<real_t>;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   SimplePAB( const BoundaryUID& boundaryUID, const FlagUID & uid, PDFField * const pdfField,
              FlagFieldT * const flagField, const real_t latticeDensity, const real_t omega, FlagUID domain, FlagUID noSlip )
      : Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfs_( pdfField ), flags_(flagField),
        latticeDensity_( latticeDensity ), omega_( omega ), tau_( real_t(1) / omega )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfs_  );
      WALBERLA_ASSERT_NOT_NULLPTR( flags_ );

      if( !flags_->flagExists( domain ) )
         domainFlag_ = flags_->registerFlag( domain );
      else
         domainFlag_ = flags_->getFlag( domain );

      if( !flags_->flagExists( noSlip ) )
         noSlipFlag_ = flags_->registerFlag( noSlip );
      else
         noSlipFlag_ = flags_->getFlag( noSlip );

      WALBERLA_ASSERT_UNEQUAL( domainFlag_, 0 );
      WALBERLA_ASSERT_UNEQUAL( noSlipFlag_, 0 );

      if( omega_ <= real_t(0) && omega_ >= real_t(2) )
      {
         WALBERLA_ABORT( "You are trying to use the simplePAB boundary condition with an omega of " << omega_ << ". "
                         "Omega has to be in the open interval (0,2)." );
      }
      else if( omega_ > real_t(0.99) && omega_ < real_t(1.01) )
      {
         WALBERLA_ABORT( "You are trying to use the simplePAB boundary condition with an omega of " << omega_ << ". "
                         "With an omega that close to 1, the pre collision PDFs can not be restored and SimplePAB will "
                         "not work!" );
      }
      else if( omega_ > real_t(0.9) && omega_ < real_t(1.1) )
      {
         WALBERLA_LOG_WARNING( "You are trying to use the simplePAB boundary condition with an omega of " << omega_ << ". "
                               "With an omega that close to 1, the pre collision PDFs can not be restored accurately and "
                               "SimplePAB will probably not work!" );
      }

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
                                                                    // current implementation of this boundary condition (SimplePAB)

      using namespace stencil;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // The cells used (Parentheses: Ginzburg notation):
      // (   x,   y,   z ): fluid cell at boundary (r_b)
      // (  nx,  ny,  nz ): boundary cell          (r_b + c_q)
      // ( nfx, nfy, nfz ): neighboring fluid cell (r_b - c_q)
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      const cell_idx_t nfx = x - cx[dir];
      const cell_idx_t nfy = y - cy[dir];
      const cell_idx_t nfz = z - cz[dir];

      //////////////////////////////////////////////////
      // Definitions of used velocities and pressures //
      //////////////////////////////////////////////////

      Vec3Real wu;  // Velocity at boundary
      Vec3Real u;   // Velocity at (x, y, z)
      Vec3Real nfu; // Velocity at (nfx, nfy, nfz)

      //////////////////////////////////////////////////
      // Get macroscopic values from PDFs             //
      //////////////////////////////////////////////////

      const real_t rho = pdfs_->getDensityAndVelocity( u, x, y, z );
                   nfu = pdfs_->getVelocity( nfx, nfy, nfz );

      //////////////////////////////////////////////////
      // Extrapolate velocity at boundary (wu)       //
      //////////////////////////////////////////////////

      if( flags_->get(nfx, nfy, nfz) & domainFlag_ )
      {
         wu = u + real_t(0.5) * (u - nfu);
      }
      else // Check whether we can get the velocity from a BC
      {
         if( flags_->get(nfx, nfy, nfz) & noSlipFlag_ )
         {
            wu = real_t(2.0) * u; // NoSlip: u at boundary is 0
         }
         else if( flags_->get(nfx, nfy, nfz) & this->mask_ )
         {
            wu = u; // PAB: u at boundary is unknown,
         }
         else
         {
            WALBERLA_ABORT( "PAB BC tries to get velocity of cell (" << nfx << ", "
                            << nfy << ", " << nfz << ") but it is is neither NoSlip nor Liquid!");
         }
      }

      ////////////////////////////////////////////////////
      // Restore pre-collision PDFs for correction term //
      ////////////////////////////////////////////////////

      real_t f1 = pdfs_->get( x, y, z, Stencil::idx[dir] );
      real_t f2 = pdfs_->get( x, y, z, Stencil::idx[inverseDir[dir]] );
      real_t f1_eq = EquilibriumDistribution<LatticeModel_T>::get(            dir , u, rho);
      real_t f2_eq = EquilibriumDistribution<LatticeModel_T>::get( inverseDir[dir], u, rho );
      real_t f1_pc = ( tau_ * f1 - f1_eq ) / ( tau_ - real_t(1) );
      real_t f2_pc = ( tau_ * f2 - f2_eq ) / ( tau_ - real_t(1) );
      real_t symDistFunc_pc = real_c( 0.5 ) * ( f1_pc + f2_pc );

      ////////////////////////////////////////////////////
      // Set boundary PDF                               //
      ////////////////////////////////////////////////////

      // result will be streamed to (x,y,z, stencil::inverseDir[d]) during sweep
      pdfs_->get(nx, ny, nz, Stencil::invDirIdx(dir)) =
         - pdfs_->get(x, y, z, Stencil::idx[dir])                                                                                    //anti-bounce-back
         + real_t(2.0) * EquilibriumDistribution<LatticeModel_T>::getSymmetricPart( dir, wu, latticeDensity_ )                       //pressure term
         + ( real_t(2.0) - omega_ ) * ( symDistFunc_pc - EquilibriumDistribution<LatticeModel_T>::getSymmetricPart( dir, u, rho ) ); //error correction
   }

protected:

   FlagUID uid_;

   PDFField   * pdfs_;
   FlagFieldT * flags_;

   flag_t noSlipFlag_;
   flag_t domainFlag_;

   real_t latticeDensity_;
   real_t omega_, tau_;

}; // class SimplePAB



} // namespace lbm
} // namespace walberla
