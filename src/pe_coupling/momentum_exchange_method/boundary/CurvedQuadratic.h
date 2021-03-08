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
//! \file CurvedQuadratic.h
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

#include <type_traits>
#include <vector>


namespace walberla {
namespace pe_coupling {


//**************************************************************************************************************************************
/*!
*   \brief Quadratic boundary handling for moving obstacles
*
*   This boundary condition implements the MR1 scheme from
*   Ginzburg et al. - "Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions" (2008)
*   References to equations and tables in the code documentation are contained in this paper.
*
*   It uses two additional cells in inverse direction of the boundary handling for a quadratic interpolation.
*   It additionally needs the pre-collision values of the PDFs for the error correction which are to be provided via the field pdfFieldPreCollision.
*   This requires the LBM sweep to be explicitly split into stream and collide to have the following form:
*   copy PDFs to extra field - collision - boundary handling - stream
*
*   In case not enough cells are available, two fall back solutions are available:
*   If only one additional cell is available, the scheme is replaced by the MGYLI scheme.
*   If no additional cell is available, the scheme is replaced by the MGYLI scheme in one-point form.
*
*   Note: The implementation could be changed if two ghost layers were used which is however omitted for efficiency reasons.
*
*/
//**************************************************************************************************************************************
template< typename LatticeModel_T, typename FlagField_T >
class CurvedQuadratic : public Boundary< typename FlagField_T::flag_t >
{
   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, lbm::collision_model::TRT_tag >::value), "Only works with TRT!" ); // to access lambda_d

   using PDFField_T = lbm::PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagField_T::flag_t;
   using BodyField_T = Field<pe::BodyID, 1>;

public:

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& )
   {
      WALBERLA_ABORT( "A CurvedQuadratic boundary cannot be created from a config file" );
      return make_shared<BoundaryConfiguration>();
   }

   inline CurvedQuadratic( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                           BodyField_T * const bodyField,  const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block,
                           PDFField_T * const pdfFieldPreCollision );

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

   PDFField_T * const pdfFieldPreCollision_;

   real_t lengthScalingFactor_;
   real_t forceScalingFactor_;

}; // class CurvedQuadratic


template< typename LatticeModel_T, typename FlagField_T>
inline CurvedQuadratic< LatticeModel_T, FlagField_T >::CurvedQuadratic( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                                                                        BodyField_T * const bodyField,  const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block,
                                                                        PDFField_T * const pdfFieldPreCollision ):
Boundary<flag_t>( boundaryUID ), uid_( uid ), domainMask_(domain), pdfField_( pdfField ), bodyField_( bodyField ), flagField_( flagField ),
blockStorage_( blockStorage ), block_( block ), pdfFieldPreCollision_ ( pdfFieldPreCollision )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField_ );
   WALBERLA_ASSERT( flagField_->isRegistered( domainMask_ )  );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfFieldPreCollision_ );

   // force scaling factor to account for different dx on this block
   const real_t dxCurrentLevel = blockStorage_.dx( blockStorage_.getLevel(block) );
   lengthScalingFactor_ = dxCurrentLevel;
   forceScalingFactor_ = lengthScalingFactor_ * lengthScalingFactor_;
}


template< typename LatticeModel_T, typename FlagField_T >
#ifndef NDEBUG
inline void CurvedQuadratic< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                            const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void CurvedQuadratic< LatticeModel_T, FlagField_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
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
   const real_t delta = lbm::intersectionRatio( body, cellCenter, direction, tolerance );

   WALBERLA_ASSERT_LESS_EQUAL( delta, real_t(1));
   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t(0));

   bool useMR1full = false;

   real_t pdf_new = real_c(0);
   real_t alpha = real_c(0);

   // get the cell indices of the cell further away from the obstacle
   const cell_idx_t xff = x - cell_idx_c( stencil::cx[ dir ] );
   const cell_idx_t yff = y - cell_idx_c( stencil::cy[ dir ] );
   const cell_idx_t zff = z - cell_idx_c( stencil::cz[ dir ] );

   auto domainWithGhostlayer = pdfField_->xyzSizeWithGhostLayer();

   // check if MR can be applied, i.e. the further-away-cell has to be a valid (i.e. fluid) cell
   if( flagField_->isPartOfMaskSet( xff, yff, zff, domainMask_ ) )
   {
      // MR1 scheme can be applied

      // get the cell indices of the cell two cells away from the obstacle
      const cell_idx_t xfff = x - cell_idx_c( 2 * stencil::cx[ dir ] );
      const cell_idx_t yfff = y - cell_idx_c( 2 * stencil::cy[ dir ] );
      const cell_idx_t zfff = z - cell_idx_c( 2 * stencil::cz[ dir ] );

      if( domainWithGhostlayer.contains( xfff, yfff, zfff ) )
      {
         if( flagField_->isPartOfMaskSet( xfff, yfff, zfff, domainMask_ ) )
         {
            useMR1full = true;
         }
      }

      if( useMR1full )
      {
         // enough cells to use the full MR1 scheme

         // coefficients from Table 3
         const real_t common = ( real_c(1) + delta ) * ( real_c(1) + delta );
         const real_t kappa0 = ( real_c(1) - real_c(2) * delta - real_c(2) * delta * delta ) / common;
         const real_t kappa_1 = delta * delta / common;
         const real_t kappaInv_1 = - kappa0;
         const real_t kappaInv_2 = - kappa_1;
         const real_t kappa1 = real_c(1) - kappa0 - kappa_1 - kappaInv_1 - kappaInv_2;

         // Eq. (4.1)
         pdf_new =   kappa1     * pdfField_->get( x   , y   , z   , Stencil_T::idx[dir] )
                   + kappa0     * pdfField_->get( xff , yff , zff , Stencil_T::idx[dir] )
                   + kappa_1    * pdfField_->get( xfff, yfff, zfff, Stencil_T::idx[dir] )
                   + kappaInv_1 * pdfField_->get( x   , y   , z   , Stencil_T::invDirIdx(dir) )
                   + kappaInv_2 * pdfField_->get( xff , yff , zff , Stencil_T::invDirIdx(dir) );

         // Table 4
         alpha = real_c(4) / common;

         // add correction term for MR1, Eq. (5.7), with pre-collision PDFs
         Vector3<real_t> velocity;
         const real_t rho = pdfFieldPreCollision_->getDensityAndEquilibriumVelocity( velocity, x, y, z );
         const real_t lambdaO = -pdfFieldPreCollision_->latticeModel().collisionModel().lambda_d();
         const real_t fqOdd = real_c(0.5) * ( pdfFieldPreCollision_->get( x, y, z, Stencil_T::idx[dir] ) - pdfFieldPreCollision_->get( x, y, z, Stencil_T::invDirIdx(dir) ) );
         const real_t nqOdd = fqOdd - lbm::EquilibriumDistribution<LatticeModel_T>::getAsymmetricPart( dir, velocity, rho );
         pdf_new += - alpha * ( real_c(0.5) * lambdaO + real_c(1) ) * nqOdd;

      }
      else
      {
         // only one additional cell available - first type special link

         // Ginzburg proposes fallback to MYLI scheme for third order, but this includes finite difference approximation in fpc (correction term)
         // here: fallback to MGYLI (formally only second order)

         // coefficients from Table 3
         const real_t common = real_c(1) + delta;
         const real_t kappa0 = ( real_c(1) - delta ) / common;
         const real_t kappaInv_1 = delta / common;
         const real_t kappa1 = real_c(1) - kappa0 - kappaInv_1;

         // Eq. (4.1)
         pdf_new =   kappa1     * pdfField_->get( x  , y  , z  , Stencil_T::idx[dir] )
                   + kappa0     * pdfField_->get( xff, yff, zff, Stencil_T::idx[dir] )
                   + kappaInv_1 * pdfField_->get( x  , y  , z  , Stencil_T::invDirIdx(dir) );

         // Table 4
         alpha = real_c(2) / common;

         // add correction term for MYLI, Eq. (5.5), with pre-collision PDFs
         Vector3<real_t> velocity;
         const real_t rho = pdfFieldPreCollision_->getDensityAndEquilibriumVelocity( velocity, x, y, z );
         const real_t lambdaO = -pdfFieldPreCollision_->latticeModel().collisionModel().lambda_d();
         const real_t fqOdd = real_c(0.5) * ( pdfFieldPreCollision_->get( x, y, z, Stencil_T::idx[dir] ) - pdfFieldPreCollision_->get( x, y, z, Stencil_T::invDirIdx(dir) ) );
         const real_t nqOdd = fqOdd - lbm::EquilibriumDistribution<LatticeModel_T>::getAsymmetricPart( dir, velocity, rho );
         const real_t mqF = lambdaO * nqOdd; // Eq. (2.17) & (2.2)
         pdf_new += - alpha * ( - real_c(0.5) ) * mqF;

      }
   }
   else
   {
      // only one fluid cell between two obstacle cells -> second type special link
      // Ginzburg proposes fallback to scheme of MLI family for third order, but this includes finite difference approximation in fpc (correction term)
      // here: fallback to MGYLI with Eq. (4.3) instead of (4.1) for one point scheme

      // coefficients from Table 3
      const real_t common = real_c(1) + delta;
      const real_t kappa0 = ( real_c(1) - delta ) / common;
      const real_t kappaInv_1 = delta / common;
      const real_t kappa1 = real_c(1) - kappa0 - kappaInv_1;

      // Eq. (4.3)
      pdf_new =   kappa1     * pdfField_->            get( x, y, z, Stencil_T::idx[dir] )
                + kappa0     * pdfFieldPreCollision_->get( x, y, z, Stencil_T::idx[dir] )
                + kappaInv_1 * pdfField_->            get( x, y, z, Stencil_T::invDirIdx(dir) );

      // Table 4
      alpha = real_c(2) / common;

      // add correction term for MYLI, Eq. (5.5), with pre-collision PDFs
      Vector3<real_t> velocity;
      const real_t rho = pdfFieldPreCollision_->getDensityAndEquilibriumVelocity( velocity, x, y, z );
      const real_t lambdaO = -pdfFieldPreCollision_->latticeModel().collisionModel().lambda_d();
      const real_t fqOdd = real_c(0.5) * ( pdfFieldPreCollision_->get( x, y, z, Stencil_T::idx[dir] ) - pdfFieldPreCollision_->get( x, y, z, Stencil_T::invDirIdx(dir) ) );
      const real_t nqOdd = fqOdd - lbm::EquilibriumDistribution<LatticeModel_T>::getAsymmetricPart( dir, velocity, rho );
      const real_t mqF = lambdaO * nqOdd; // Eq. (2.17) & (2.2)
      pdf_new += - alpha * ( - real_c(0.5) ) * mqF;

   }

   // get coordinates of fluid cell center
   real_t cx, cy, cz;
   blockStorage_.getBlockLocalCellCenter( block_, Cell(x,y,z), cx, cy, cz );

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
