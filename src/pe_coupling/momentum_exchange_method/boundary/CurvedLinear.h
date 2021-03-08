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
//! \file CurvedLinear.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/all.h"

#include "field/FlagField.h"

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/ForceModel.h"

#include "stencil/Directions.h"

#include "pe_coupling/geometry/PeIntersectionRatio.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/Types.h"
#include <vector>


namespace walberla {
namespace pe_coupling {

//**************************************************************************************************************************************
/*!
*   \brief Linear boundary handling for moving obstacles
*
*   This boundary condition implements the CLI scheme from
*   Ginzburg et al. - "Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions" (2008)
*   References to equations and tables in the code documentation are contained in this paper.
*
*   It uses one additional cell in inverse direction of the boundary handling for a linear interpolation.
*
*   In case not enough cells are available, one fall back solution is available:
*   If no additional cell is available, the scheme is replaced by the bounce back scheme (see SimpleBB.h).
*
*/
//**************************************************************************************************************************************
template< typename LatticeModel_T, typename FlagField_T >
class CurvedLinear : public Boundary< typename FlagField_T::flag_t >
{
   using PDFField_T = lbm::PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagField_T::flag_t;
   using BodyField_T = Field<pe::BodyID, 1>;

public:

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& )
   {
      WALBERLA_ABORT( "A CurvedLinear boundary cannot be created from a config file" );
      return make_shared<BoundaryConfiguration>();
   }

   inline CurvedLinear( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                        BodyField_T * const bodyField,  const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block );

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

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

private:

   const FlagUID uid_;

   flag_t domainMask_;

   PDFField_T  * const  pdfField_;
   BodyField_T * const bodyField_;
   const FlagField_T * const flagField_;

   const StructuredBlockStorage & blockStorage_;
   const IBlock & block_;

   real_t lengthScalingFactor_;
   real_t forceScalingFactor_;

}; // class CurvedLinear


template< typename LatticeModel_T, typename FlagField_T >
inline CurvedLinear< LatticeModel_T, FlagField_T >::CurvedLinear( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                                                                  BodyField_T * const bodyField,  const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block ):
Boundary<flag_t>( boundaryUID ), uid_( uid ), domainMask_(domain), pdfField_( pdfField ), bodyField_( bodyField ), flagField_( flagField ), blockStorage_( blockStorage ), block_( block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField_ );
   WALBERLA_ASSERT( flagField_->isRegistered( domainMask_ )  );

   // force scaling factor to account for different dx on this block
   const real_t dxCurrentLevel = blockStorage_.dx( blockStorage_.getLevel(block) );
   lengthScalingFactor_ = dxCurrentLevel;
   forceScalingFactor_ = lengthScalingFactor_ * lengthScalingFactor_;
}


template< typename LatticeModel_T, typename FlagField_T >
#ifndef NDEBUG
inline void CurvedLinear< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void CurvedLinear< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL  ( nx, x + cell_idx_t( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL  ( ny, y + cell_idx_t( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL  ( nz, z + cell_idx_t( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL  ( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                               // current implementation of this boundary condition
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField_->get(nx,ny,nz) );

   // determine distance to real boundary, i.e. delta value
   // cell center of the near-boundary fluid cell
   Cell nearBoundaryCell(x,y,z);
   Vector3< real_t > cellCenter = blockStorage_.getBlockLocalCellCenter(block_, nearBoundaryCell);

   // direction of the ray (from the fluid cell center to the boundary cell)
   Vector3< real_t > direction( lengthScalingFactor_ * real_c( stencil::cx[ dir ] ),
                                lengthScalingFactor_ * real_c( stencil::cy[ dir ] ),
                                lengthScalingFactor_ * real_c( stencil::cz[ dir ] ) );

   // (if applicable) line search accuracy
   const real_t tolerance = real_c( 1e-4 ) * lengthScalingFactor_;

   //get Body
   pe::RigidBody & body = *(bodyField_->get( nx, ny, nz ));

   // delta value obtained by ray - body intersection
   // depending on the implementation for the specific body, either an analytical formula (e.g. for the sphere) or a line search algorithm is used
   real_t delta = lbm::intersectionRatio( body, cellCenter, direction, tolerance );

   WALBERLA_ASSERT_LESS_EQUAL( delta, real_t(1));
   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t(0));

   real_t pdf_new = real_t(0);
   real_t alpha = real_t(0);

   // get the cell indices of the cell further away from the obstacle
   const cell_idx_t xff = x - cell_idx_c( stencil::cx[ dir ] );
   const cell_idx_t yff = y - cell_idx_c( stencil::cy[ dir ] );
   const cell_idx_t zff = z - cell_idx_c( stencil::cz[ dir ] );

   // check if CLI can be applied, i.e. the further-away-cell has to be a valid (i.e. fluid) cell
   if( flagField_->isPartOfMaskSet( xff, yff, zff, domainMask_ ) )
   {
      // apply CLI scheme

      // coefficients from Table 3
      const real_t kappa0 = ( real_c(1) - real_c(2) * delta ) / ( real_c(1) + real_c(2) * delta );

      // Eq. (4.1)
      pdf_new = /*kappa1*/ pdfField_->get( x  , y  , z  , Stencil_T::idx[dir] )
                + kappa0 * pdfField_->get( xff, yff, zff, Stencil_T::idx[dir] )
                - kappa0 * pdfField_->get( x  , y  , z  , Stencil_T::invDirIdx(dir) );

      // Table 4
      alpha = real_c(4) / ( real_c(1) + real_c(2) * delta);

   }
   else
   {
      // fall back to Bounce back

      delta = real_c(0.5);
      pdf_new = pdfField_->get( x, y, z, Stencil_T::idx[dir] );
      alpha = real_c(2);
   }

   // assumed boundary position
   const Vector3<real_t> boundaryPosition( cellCenter + delta * direction );

   // obtain the velocity at the obstacle boundary
   auto boundaryVelocity = body.velFromWF( boundaryPosition );

   // include effect of boundary velocity
   if( LatticeModel_T::compressible )
   {
      const auto density  = pdfField_->getDensity(x,y,z);
      pdf_new -= real_c(3.0) * alpha * density * LatticeModel_T::w[ Stencil_T::idx[dir] ] *
                 ( real_c( stencil::cx[ dir ] ) * boundaryVelocity[0] +
                   real_c( stencil::cy[ dir ] ) * boundaryVelocity[1] +
                   real_c( stencil::cz[ dir ] ) * boundaryVelocity[2] );
   }
   else
   {
      pdf_new -= real_c(3.0) * alpha * LatticeModel_T::w[ Stencil_T::idx[dir] ] *
                 ( real_c( stencil::cx[ dir ] ) * boundaryVelocity[0] +
                   real_c( stencil::cy[ dir ] ) * boundaryVelocity[1] +
                   real_c( stencil::cz[ dir ] ) * boundaryVelocity[2] );
   }

   // carry out the boundary handling
   pdfField_->get( nx, ny, nz, Stencil_T::invDirIdx(dir) ) = pdf_new;

   // check if fluid cell is not inside ghost layer ( = is in inner part), since then no force is allowed to be added
   if( !pdfField_->isInInnerPart( nearBoundaryCell ) )
   {
      return;
   }

   // calculate the force on the obstacle
   // original work (MEM):      Ladd - Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation (1994)
   // improved version (GIMEM): Wen et al. - Galilean invariant fluid-solid interfacial dynamics in lattice Boltzmann simulations (2014)

   const real_t pdf_old = pdfField_->get( x, y, z, Stencil_T::idx[dir] );

   // MEM: F = pdf_old + pdf_new - common
   const real_t forceMEM = pdf_old + pdf_new;

   // correction from Wen
   const real_t correction = pdf_old - pdf_new;

   // force consists of the MEM part and the galilean invariance correction including the boundary velocity
   Vector3<real_t> force( real_c( stencil::cx[dir] ) * forceMEM - correction * boundaryVelocity[0],
                          real_c( stencil::cy[dir] ) * forceMEM - correction * boundaryVelocity[1],
                          real_c( stencil::cz[dir] ) * forceMEM - correction * boundaryVelocity[2] );

   force *= forceScalingFactor_;

   // add the force onto the body at the obstacle boundary
   body.addForceAtPos( force, boundaryPosition );

   /*
   WALBERLA_LOG_DETAIL_SECTION() {
      std::stringstream ss;
      ss << "MOBoundary in cell <" << x << ", " << y << ", " << z << "> in dir <"
         << stencil::cx[dir] << ", " << stencil::cy[dir] << ", " << stencil::cz[dir]
         << ">:  on body with id " << bodyField_->get(nx,ny,nz)->getSystemID()
         << ", applying force " << force << " at position " << boundaryPosition
         << ", with pdf_old=" << pdf_old << ", pdf_new=" << pdf_new;
      WALBERLA_LOG_DETAIL( ss.str() );
   }
   */
}

} // namespace pe_coupling
} // namespace walberla
