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
//! \file Reconstructor.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/Field.h"
#include "lbm/field/PdfField.h"
#include "ExtrapolationDirectionFinder.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/Types.h"

#include <type_traits>


namespace walberla {
namespace pe_coupling {



//**************************************************************************************************************************************
/*!
*   \brief Classes to be used together with the PDFReconstruction class to reconstruct PDFs
*
*   Each reconstructor must exactly implement the member function
*     void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );
*   that reconstructs all PDFs in a specific cell with indices x,y,z on a given block.
*
*   Different variants are available:
*
*   - EquilibriumReconstructor:
*     Determines a local average density and sets the PDFs to their equilibrium values based on this density and the body's velocity
*
*   - EquilibriumAndNonEquilibriumReconstructor:
*     First reconstructs the equilibrium values with the help of the EquilibriumReconstructor.
*     Then extrapolates the non-equilibrium part from a neighboring cell and adds it to the reconstructed PDFs for better accuracy.
*     Defaults to EquilibriumReconstructor if no suitable cell for extrapolation is available.
*
*   - ExtrapolationReconstructor:
*     Extrapolates the PDFs of three or two neighboring cells that lie in extrapolation direction.
*     Optionally, a special treatment after the extrapolation can be used that sets certain moments directly (only D3Q19).
*     Defaults to EquilibriumAndNonEquilibriumReconstructor if not enough cells in this direction are available.
*
*   EquilibriumAndNonEquilibriumReconstructor and ExtrapolationReconstructor need an extrapolation direction
*   which is provided by the ExtrapolationDirectionFinder that is to be chosen (see ExtrapolationDirectionFinder.h).
*
*/
//**************************************************************************************************************************************


template< typename LatticeModel_T, typename BoundaryHandling_T >
class EquilibriumReconstructor
{
public:

   typedef lbm::PdfField< LatticeModel_T > PdfField_T;
   typedef Field< pe::BodyID, 1 >          BodyField_T;

   EquilibriumReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                             const BlockDataID & bodyFieldID )
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ), bodyFieldID_( bodyFieldID )
   {}

   void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );

private:

   void setLocalEquilibrium( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );

   real_t getLocalAverageDensity( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField ) const;

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyFieldID_;

};

template< typename LatticeModel_T, typename BoundaryHandling_T >
void EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T >
::operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   setLocalEquilibrium( x, y, z, block, pdfField);
}

template< typename LatticeModel_T, typename BoundaryHandling_T >
void EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T >
::setLocalEquilibrium( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BodyField_T * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );
   WALBERLA_ASSERT_EQUAL( pdfField->xyzSize(), bodyField->xyzSize() );

   const real_t averageDensity = getLocalAverageDensity( x, y, z, block, pdfField );

   real_t cx, cy, cz;
   blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

   WALBERLA_ASSERT_NOT_NULLPTR( (*bodyField)(x,y,z) );
   const auto velocity = (*bodyField)(x,y,z)->velFromWF(cx,cy,cz);

   pdfField->setToEquilibrium( x, y, z, Vector3< real_t >( velocity[0], velocity[1], velocity[2] ), averageDensity );
}

template< typename LatticeModel_T, typename BoundaryHandling_T >
real_t EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T >
::getLocalAverageDensity( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

   CellInterval localDomain = pdfField->xyzSize();

   uint_t nAverage = uint_t(0);
   real_t averageDensity = real_t(0);
   for( auto neighborDir = stencil::D3Q27::beginNoCenter(); neighborDir != stencil::D3Q27::end(); ++neighborDir )
   {
      Cell neighbor( x + neighborDir.cx(), y + neighborDir.cy(), z + neighborDir.cz() );

      // check if neighbor cell is inside domain
      if( !( localDomain.contains( neighbor ) ) )
         continue;

      // check if neighbor cell is a fluid cell
      if( boundaryHandling->isDomain( neighbor ) )
      {
         // obtain density and add it up
         averageDensity += pdfField->getDensity( neighbor );
         ++nAverage;
      }
   }

   return ( nAverage > uint_t( 0 ) ) ? averageDensity / real_c( nAverage ) : real_t(1.0);
}


template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
class EquilibriumAndNonEquilibriumReconstructor
{
public:

   typedef lbm::PdfField< LatticeModel_T > PdfField_T;
   typedef Field< pe::BodyID, 1 >          BodyField_T;

   EquilibriumAndNonEquilibriumReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                                              const BlockDataID & bodyFieldID,
                                              const ExtrapolationDirectionFinder_T & extrapolationDirectionFinder )
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ),
     bodyFieldID_( bodyFieldID ),
     extrapolationDirectionFinder_( extrapolationDirectionFinder ),
     equilibriumReconstructor_( EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T >( blockStorage, boundaryHandlingID, bodyFieldID ) )
   { }

   void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );
   void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField, const Vector3<cell_idx_t> & extrapolationDirection );

private:

   void setNonEquilibrium( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField, const Vector3<cell_idx_t> & extrapolationDirection );

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyFieldID_;

   ExtrapolationDirectionFinder_T extrapolationDirectionFinder_;
   EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T > equilibriumReconstructor_;

};

template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   Vector3<cell_idx_t> extrapolationDirection(0);
   extrapolationDirectionFinder_.getDirection( x, y, z, block, extrapolationDirection );

   (*this)(x, y, z, block, pdfField, extrapolationDirection);
}

template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField, const Vector3<cell_idx_t> & extrapolationDirection  )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   equilibriumReconstructor_( x, y, z, block, pdfField );

   if( extrapolationDirection != cell_idx_t(0) )
   {
      setNonEquilibrium( x, y, z, block, pdfField, extrapolationDirection);
   }
}

template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::setNonEquilibrium( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField, const Vector3<cell_idx_t> & extrapolationDirection )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

   CellInterval localDomain = pdfField->xyzSize();
   Cell normalNeighbor( x + extrapolationDirection[0], y + extrapolationDirection[1], z + extrapolationDirection[2] );
   if( localDomain.contains( normalNeighbor ) )
   {
      if( boundaryHandling->isDomain( normalNeighbor ) )
      {
         real_t neighborDensity = pdfField->getDensity( normalNeighbor );
         Vector3< real_t > neighborVelocity = pdfField->getVelocity( normalNeighbor );
         std::vector< real_t > equilibriumDistributions = lbm::EquilibriumDistribution< LatticeModel_T >::get( neighborVelocity, neighborDensity );

         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
         {
            pdfField->get( x, y, z, d.toIdx() ) += 
                     ( pdfField->get( x + extrapolationDirection[0], y + extrapolationDirection[1], z + extrapolationDirection[2], d.toIdx() ) - equilibriumDistributions[ d.toIdx() ] );
         }
      }
   }
}


template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
class ExtrapolationReconstructor
{
public:

   typedef lbm::PdfField< LatticeModel_T > PdfField_T;
   typedef Field< pe::BodyID, 1 >          BodyField_T;

   ExtrapolationReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                               const BlockDataID & bodyFieldID,
                               const ExtrapolationDirectionFinder_T & extrapolationDirectionFinder,
                               const bool & enforceNoSlipConstraintAfterExtrapolation = false )
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ),
     bodyFieldID_( bodyFieldID ),
     enforceNoSlipConstraintAfterExtrapolation_( enforceNoSlipConstraintAfterExtrapolation ),
     extrapolationDirectionFinder_( extrapolationDirectionFinder ),
     alternativeReconstructor_( EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
      ( blockStorage, boundaryHandlingID, bodyFieldID, extrapolationDirectionFinder ) )
   {
      if( enforceNoSlipConstraintAfterExtrapolation_ ) {
         WALBERLA_CHECK((std::is_same<typename LatticeModel_T::Stencil, stencil::D3Q19>::value),
                        "Enforcing no-slip constraint after extrapolation currently only works with D3Q19 stencil!");
      }
   }

   void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );

private:

   uint_t getNumberOfExtrapolationCells( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField,
                                         const Vector3<cell_idx_t> & extrapolationDirection ) const;

   void extrapolatePDFs( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField,
                         const Vector3<cell_idx_t> & extrapolationDirection, const uint_t & numberOfCellsForExtrapolation);

   void enforceNoSlipConstraint( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField );

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyFieldID_;

   const bool enforceNoSlipConstraintAfterExtrapolation_;

   ExtrapolationDirectionFinder_T extrapolationDirectionFinder_;
   EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T > alternativeReconstructor_;

};


template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void ExtrapolationReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   Vector3<cell_idx_t> extrapolationDirection(0);
   extrapolationDirectionFinder_.getDirection( x, y, z, block, extrapolationDirection );

   uint_t numberOfCellsForExtrapolation = getNumberOfExtrapolationCells( x, y, z, block, pdfField, extrapolationDirection );

   if( numberOfCellsForExtrapolation < uint_t(2) )
   {
      alternativeReconstructor_( x, y, z, block, pdfField, extrapolationDirection );
   }
   else
   {
      extrapolatePDFs( x, y, z, block, pdfField, extrapolationDirection, numberOfCellsForExtrapolation );
      if( enforceNoSlipConstraintAfterExtrapolation_ )
      {
         enforceNoSlipConstraint( x, y, z, block, pdfField );
      }
   }
}

template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
uint_t ExtrapolationReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::getNumberOfExtrapolationCells( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField,
                                 const Vector3<cell_idx_t> & extrapolationDirection ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

   if( extrapolationDirection == cell_idx_t(0) ) return uint_t(0);

   CellInterval localDomain = pdfField->xyzSize();

   uint_t desiredCellsInExtrapolationDirection = uint_t(3);

   for( uint_t numCells = uint_t(1); numCells <= desiredCellsInExtrapolationDirection; ++numCells )
   {
      Cell checkCell( x + cell_idx_c(numCells) * extrapolationDirection[0], y + cell_idx_c(numCells) * extrapolationDirection[1], z + cell_idx_c(numCells) * extrapolationDirection[2] );

      // check if cell is inside domain & fluid
      if( !localDomain.contains( checkCell ) )
         return numCells - uint_t(1);

      if( !boundaryHandling->isDomain( checkCell ) )
         return numCells - uint_t(1);
   }
   return desiredCellsInExtrapolationDirection;
}


template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void ExtrapolationReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::extrapolatePDFs( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const /*block*/, PdfField_T * const pdfField,
                   const Vector3<cell_idx_t> & extrapolationDirection, const uint_t & numberOfCellsForExtrapolation)
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

   if( numberOfCellsForExtrapolation == uint_t(3) )
   {
      // quadratic normal extrapolation
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         pdfField->get( x, y, z, d.toIdx() ) = real_t(3) * pdfField->get( x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], d.toIdx() )
                                             - real_t(3) * pdfField->get( x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], d.toIdx() )
                                             +             pdfField->get( x + 3*extrapolationDirection[0], y + 3*extrapolationDirection[1], z + 3*extrapolationDirection[2], d.toIdx() );
      }
   } else { // numberOfCellsForExtrapolation == 2
      // linear normal extrapolation
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         pdfField->get( x, y, z, d.toIdx() ) = real_t(2) * pdfField->get( x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], d.toIdx() )
                                             -             pdfField->get( x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], d.toIdx() );
      }
   }
}

template< typename LatticeModel_T, typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
void ExtrapolationReconstructor< LatticeModel_T, BoundaryHandling_T, ExtrapolationDirectionFinder_T >
::enforceNoSlipConstraint( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, PdfField_T * const pdfField )
{
   //NOTE: this currently works only for D3Q19 stencils!

   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BodyField_T * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );
   WALBERLA_ASSERT_EQUAL( pdfField->xyzSize(), bodyField->xyzSize() );

   real_t cx, cy, cz;
   blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

   WALBERLA_ASSERT_NOT_NULLPTR( (*bodyField)(x,y,z) );
   Vector3<real_t> bodyVelocity = (*bodyField)(x,y,z)->velFromWF(cx,cy,cz);
   WALBERLA_ASSERT( !math::isnan(bodyVelocity) );

   // transforms to moment space (see MRT collision model) to set the body's velocity in cell without affecting other moments
   const real_t _1_2  = real_t(1) / real_t(2);
   const real_t _1_3  = real_t(1) / real_t(3);
   const real_t _1_4  = real_t(1) / real_t(4);
   const real_t _1_6  = real_t(1) / real_t(6);
   const real_t _1_8  = real_t(1) / real_t(8);
   const real_t _1_12 = real_t(1) / real_t(12);
   const real_t _1_16 = real_t(1) / real_t(16);
   const real_t _1_18 = real_t(1) / real_t(18);
   const real_t _1_24 = real_t(1) / real_t(24);
   const real_t _1_36 = real_t(1) / real_t(36);
   const real_t _1_48 = real_t(1) / real_t(48);
   const real_t _1_72 = real_t(1) / real_t(72);

   // restriction to D3Q19!
   const real_t vC  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::C]  );
   const real_t vN  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::N]  );
   const real_t vS  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::S]  );
   const real_t vW  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::W]  );
   const real_t vE  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::E]  );
   const real_t vT  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::T]  );
   const real_t vB  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::B]  );
   const real_t vNW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NW] );
   const real_t vNE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NE] );
   const real_t vSW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SW] );
   const real_t vSE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SE] );
   const real_t vTN = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TN] );
   const real_t vTS = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TS] );
   const real_t vTW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TW] );
   const real_t vTE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TE] );
   const real_t vBN = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BN] );
   const real_t vBS = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BS] );
   const real_t vBW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BW] );
   const real_t vBE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BE] );

   // transform to moment space and change momentum to body's velocity ( * rho_0 )
   const real_t m0  = vC + vN + vS + vW + vE + vT + vB + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
   const real_t m1  = -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
   const real_t m2  = vC - real_t(2) * ( vN + vS + vW + vE + vT + vB ) + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
   const real_t m3  = bodyVelocity[0];
   const real_t m4  = real_t(2) * vW - real_t(2) * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE;
   const real_t m5  = bodyVelocity[1];
   const real_t m6  = real_t(-2) * vN + real_t(2) * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS;
   const real_t m7  = bodyVelocity[2];
   const real_t m8  = real_t(-2) * vT + real_t(2) * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE;
   const real_t m9  = -vN - vS + real_t(2) * vW + real_t(2) * vE - vT - vB + vNW + vNE + vSW + vSE - real_t(2) * vTN
                      - real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE;
   const real_t m10 = vN + vS - real_t(2) * vW - real_t(2) * vE + vT + vB + vNW + vNE + vSW + vSE - real_t(2) * vTN
                      - real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE;
   const real_t m11 = vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE;
   const real_t m12 = -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE;
   const real_t m13 = -vNW + vNE + vSW - vSE;
   const real_t m14 =  vTN - vTS - vBN + vBS;
   const real_t m15 = -vTW + vTE + vBW - vBE;
   const real_t m16 = -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE;
   const real_t m17 = -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS;
   const real_t m18 = -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE;

   // transform back to velocity space
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::C] )  = _1_3  * m0  - _1_2  * m1  + _1_6  * m2;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::N] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m5  - _1_6  * m6  - _1_24 * m9  +
                                                                         _1_24 * m10 + _1_8  * m11 - _1_8  * m12;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::S] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m5  + _1_6  * m6  - _1_24 * m9  +
                                                                         _1_24 * m10 + _1_8  * m11 - _1_8  * m12;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::W] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m3  + _1_6  * m4  + _1_12 * m9  - _1_12 * m10;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::E] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m3  - _1_6  * m4  + _1_12 * m9  - _1_12 * m10;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::T] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m7  - _1_6  * m8  - _1_24 * m9  +
                                                                         _1_24 * m10 - _1_8  * m11 + _1_8  * m12;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::B] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m7  + _1_6  * m8  - _1_24 * m9  +
                                                                         _1_24 * m10 - _1_8  * m11 + _1_8  * m12;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  +
                                                                         _1_12 * m5  + _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                         _1_16 * m12 - _1_4  * m13 - _1_8  * m16 - _1_8  * m17;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  +
                                                                         _1_12 * m5  + _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                         _1_16 * m12 + _1_4  * m13 + _1_8  * m16 - _1_8  * m17;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  -
                                                                         _1_12 * m5  - _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                         _1_16 * m12 + _1_4  * m13 - _1_8  * m16 + _1_8  * m17;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  -
                                                                         _1_12 * m5  - _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                         _1_16 * m12 - _1_4  * m13 + _1_8  * m16 + _1_8  * m17;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TN] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m5  + _1_24 * m6  +
                                                                         _1_12 * m7  + _1_24 * m8  - _1_24 * m9  - _1_24 * m10 + _1_4  * m14 +
                                                                         _1_8  * m17 - _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TS] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m5  - _1_24 * m6  +
                                                                         _1_12 * m7  + _1_24 * m8  - _1_24 * m9  - _1_24 * m10 - _1_4  * m14 -
                                                                         _1_8  * m17 - _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  +
                                                                         _1_12 * m7  + _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                         _1_16 * m12 - _1_4  * m15 + _1_8  * m16 + _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  +
                                                                         _1_12 * m7  + _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                         _1_16 * m12 + _1_4  * m15 - _1_8  * m16 + _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BN] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m5  + _1_24 * m6  -
                                                                         _1_12 * m7  - _1_24 * m8  - _1_24 * m9  - _1_24 * m10 - _1_4  * m14 +
                                                                         _1_8  * m17 + _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BS] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m5  - _1_24 * m6  -
                                                                         _1_12 * m7  - _1_24 * m8  - _1_24 * m9  - _1_24 * m10 + _1_4  * m14 -
                                                                         _1_8  * m17 + _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  -
                                                                         _1_12 * m7  - _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                         _1_16 * m12 + _1_4  * m15 + _1_8  * m16 - _1_8  * m18;
   pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  -
                                                                         _1_12 * m7  - _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                         _1_16 * m12 - _1_4  * m15 - _1_8  * m16 - _1_8  * m18;
}


} // namespace pe_coupling
} // namespace walberla
