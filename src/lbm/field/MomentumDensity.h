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
//! \file MomentumDensity.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <type_traits>


// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

namespace internal {

template< typename LatticeModel_T, typename FieldPtrOrIterator >
void getMomentumDensity( Vector3< real_t > & momentumDensity, const FieldPtrOrIterator & it )
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false),
                  "If you have the D3Q19 stencil, you should use the optimized function!" );

   auto d = LatticeModel_T::Stencil::begin();

   const auto & firstPdf = it[ d.toIdx() ];

   momentumDensity[0] = firstPdf * real_c( d.cx() );
   momentumDensity[1] = firstPdf * real_c( d.cy() );
   momentumDensity[2] = firstPdf * real_c( d.cz() );

   ++d;

   while( d != LatticeModel_T::Stencil::end() )
   {
      const auto & pdf = it[ d.toIdx() ];

      momentumDensity[0] += pdf * real_c( d.cx() );
      momentumDensity[1] += pdf * real_c( d.cy() );
      momentumDensity[2] += pdf * real_c( d.cz() );

      ++d;
   }
}

template< typename LatticeModel_T, typename PdfField_T >
void getMomentumDensity( Vector3< real_t > & momentumDensity, const PdfField_T & pdf,
                         const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false),
                  "If you have the D3Q19 stencil, you should use the optimized function!" );

   const auto & xyz0 = pdf(x,y,z,0);

   auto d = LatticeModel_T::Stencil::begin();

   const auto & firstPdf = pdf.getF( &xyz0, d.toIdx() );

   momentumDensity[0] = firstPdf * real_c( d.cx() );
   momentumDensity[1] = firstPdf * real_c( d.cy() );
   momentumDensity[2] = firstPdf * real_c( d.cz() );

   ++d;

   while( d != LatticeModel_T::Stencil::end() )
   {
      const auto & pdfValue = pdf.getF( &xyz0, d.toIdx() );

      momentumDensity[0] += pdfValue * real_c( d.cx() );
      momentumDensity[1] += pdfValue * real_c( d.cy() );
      momentumDensity[2] += pdfValue * real_c( d.cz() );

      ++d;
   }
}

template< typename LatticeModel_T, typename FieldPtrOrIterator >
void getMomentumDensityD3Q19( Vector3< real_t > & momentumDensity, const FieldPtrOrIterator & it )
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "This function works with D3Q19 stencils only!" );

   using namespace stencil;
   using Stencil = typename LatticeModel_T::Stencil;

   momentumDensity[0] = it[ Stencil::idx[E] ] + it[ Stencil::idx[NE] ] + it[ Stencil::idx[SE] ] + it[ Stencil::idx[TE] ] + it[ Stencil::idx[BE] ] -
                        it[ Stencil::idx[W] ] - it[ Stencil::idx[NW] ] - it[ Stencil::idx[SW] ] - it[ Stencil::idx[TW] ] - it[ Stencil::idx[BW] ];
   momentumDensity[1] = it[ Stencil::idx[N] ] + it[ Stencil::idx[NW] ] + it[ Stencil::idx[TN] ] + it[ Stencil::idx[BN] ] + it[ Stencil::idx[NE] ] -
                        it[ Stencil::idx[S] ] - it[ Stencil::idx[SW] ] - it[ Stencil::idx[SE] ] - it[ Stencil::idx[TS] ] - it[ Stencil::idx[BS] ];
   momentumDensity[2] = it[ Stencil::idx[T] ] + it[ Stencil::idx[TS] ] + it[ Stencil::idx[TW] ] + it[ Stencil::idx[TN] ] + it[ Stencil::idx[TE] ] -
                        it[ Stencil::idx[B] ] - it[ Stencil::idx[BN] ] - it[ Stencil::idx[BS] ] - it[ Stencil::idx[BW] ] - it[ Stencil::idx[BE] ];
}

template< typename LatticeModel_T, typename PdfField_T >
void getMomentumDensityD3Q19( Vector3< real_t > & momentumDensity, const PdfField_T & pdf,
                              const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "This function works with D3Q19 stencils only!" );

   using namespace stencil;
   using Stencil = typename LatticeModel_T::Stencil;

   const auto & xyz0 = pdf(x,y,z,0);

   momentumDensity[0] = pdf.getF( &xyz0, Stencil::idx[E]  ) + pdf.getF( &xyz0, Stencil::idx[NE] ) + pdf.getF( &xyz0, Stencil::idx[SE] ) +
                        pdf.getF( &xyz0, Stencil::idx[TE] ) + pdf.getF( &xyz0, Stencil::idx[BE] ) - pdf.getF( &xyz0, Stencil::idx[W]  ) -
                        pdf.getF( &xyz0, Stencil::idx[NW] ) - pdf.getF( &xyz0, Stencil::idx[SW] ) - pdf.getF( &xyz0, Stencil::idx[TW] ) -
                        pdf.getF( &xyz0, Stencil::idx[BW] );
   momentumDensity[1] = pdf.getF( &xyz0, Stencil::idx[N]  ) + pdf.getF( &xyz0, Stencil::idx[NW] ) + pdf.getF( &xyz0, Stencil::idx[TN] ) +
                        pdf.getF( &xyz0, Stencil::idx[BN] ) + pdf.getF( &xyz0, Stencil::idx[NE] ) - pdf.getF( &xyz0, Stencil::idx[S]  ) -
                        pdf.getF( &xyz0, Stencil::idx[SW] ) - pdf.getF( &xyz0, Stencil::idx[SE] ) - pdf.getF( &xyz0, Stencil::idx[TS] ) -
                        pdf.getF( &xyz0, Stencil::idx[BS] );
   momentumDensity[2] = pdf.getF( &xyz0, Stencil::idx[T]  ) + pdf.getF( &xyz0, Stencil::idx[TS] ) + pdf.getF( &xyz0, Stencil::idx[TW] ) +
                        pdf.getF( &xyz0, Stencil::idx[TN] ) + pdf.getF( &xyz0, Stencil::idx[TE] ) - pdf.getF( &xyz0, Stencil::idx[B]  ) -
                        pdf.getF( &xyz0, Stencil::idx[BN] ) - pdf.getF( &xyz0, Stencil::idx[BS] ) - pdf.getF( &xyz0, Stencil::idx[BW] ) -
                        pdf.getF( &xyz0, Stencil::idx[BE] );
}

/////////////////////////////////////////////
// force correction - macroscopic velocity //
/////////////////////////////////////////////

template< typename LatticeModel_T, class Enable = void >
struct MacroscopicForceCorrection
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::internal::MacroscopicForceCorrection' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T >
struct MacroscopicForceCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   static void apply( const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity();
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator &, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      apply( latticeModel, momentumDensity );
   }

   static void apply( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      apply( latticeModel, momentumDensity );
   }
};

template< typename LatticeModel_T >
struct MacroscopicForceCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   /*
   static void apply( const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity();
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }
   */

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity(it.x(), it.y(), it.z());
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }

   static void apply( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity(x,y,z);
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }
};

template< typename LatticeModel_T >
struct MacroscopicForceCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::shiftMacVel >::type >
{
   static void apply( const LatticeModel_T &, Vector3< real_t > & ) {}

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator &, const LatticeModel_T &, Vector3< real_t > & ) {}

   static void apply( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T &, Vector3< real_t > & ) {}
};

/////////////////////////////////////////////
// force correction - equilibrium velocity //
/////////////////////////////////////////////

template< typename LatticeModel_T, class Enable = void >
struct EquilibriumForceCorrection
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::internal::EquilibriumForceCorrection' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T >
struct EquilibriumForceCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::ForceModel::constant &&
	                                                                         LatticeModel_T::ForceModel::shiftEquVel
	                                                                         >::type >
{
   static void apply( const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity();
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator &, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      apply( latticeModel, momentumDensity );
   }

   static void apply( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      apply( latticeModel, momentumDensity );
   }
};

template< typename LatticeModel_T >
struct EquilibriumForceCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   /*
   static void apply( const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity();
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }
   */

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity(it.x(), it.y(), it.z());
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }

   static void apply( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel, Vector3< real_t > & momentumDensity )
   {
      const auto & force = latticeModel.forceModel().forceDensity(x,y,z);
      const real_t dt_2 = real_t(0.5);

      momentumDensity[0] += dt_2 * force[0];
      momentumDensity[1] += dt_2 * force[1];
      momentumDensity[2] += dt_2 * force[2];
   }
};

template< typename LatticeModel_T >
struct EquilibriumForceCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::shiftEquVel >::type >
{
   static void apply( const LatticeModel_T &, Vector3< real_t > & ) {}

   template< typename FieldPtrOrIterator >
   static void apply( FieldPtrOrIterator &, const LatticeModel_T &, Vector3< real_t > & ) {}

   static void apply( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T &, Vector3< real_t > & ) {}
};

} // namespace internal



template< typename LatticeModel_T, class Enable = void >
struct MomentumDensity
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::MomentumDensity' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T >
struct MomentumDensity< LatticeModel_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "For D3Q19 there is an optimized version!" );

   template< typename FieldPtrOrIterator >
   static void getEquilibrium( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
   {
      internal::getMomentumDensity< LatticeModel_T, FieldPtrOrIterator >( momentumDensity, it );
      internal::EquilibriumForceCorrection< LatticeModel_T >::apply( it, latticeModel, momentumDensity );
   }

   template< typename PdfField_T >
   static void getEquilibrium( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                               const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::getMomentumDensity< LatticeModel_T, PdfField_T >( momentumDensity, pdf, x, y, z );
      internal::EquilibriumForceCorrection< LatticeModel_T >::apply( x, y, z, latticeModel, momentumDensity );
   }

   template< typename FieldPtrOrIterator >
   static void get( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
   {
      internal::getMomentumDensity< LatticeModel_T, FieldPtrOrIterator >( momentumDensity, it );
      internal::MacroscopicForceCorrection< LatticeModel_T >::apply( it, latticeModel, momentumDensity );
   }

   template< typename PdfField_T >
   static void get( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::getMomentumDensity< LatticeModel_T, PdfField_T >( momentumDensity, pdf, x, y, z );
      internal::MacroscopicForceCorrection< LatticeModel_T >::apply( x, y, z, latticeModel, momentumDensity );
   }
};

template< typename LatticeModel_T >
struct MomentumDensity< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );

   template< typename FieldPtrOrIterator >
   static void getEquilibrium( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
   {
      internal::getMomentumDensityD3Q19< LatticeModel_T, FieldPtrOrIterator >( momentumDensity, it );
      internal::EquilibriumForceCorrection< LatticeModel_T >::apply( it, latticeModel, momentumDensity );
   }

   template< typename PdfField_T >
   static void getEquilibrium( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                               const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::getMomentumDensityD3Q19< LatticeModel_T, PdfField_T >( momentumDensity, pdf, x, y, z );
      internal::EquilibriumForceCorrection< LatticeModel_T >::apply( x, y, z, latticeModel, momentumDensity );
   }

   template< typename FieldPtrOrIterator >
   static void get( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
   {
      internal::getMomentumDensityD3Q19< LatticeModel_T, FieldPtrOrIterator >( momentumDensity, it );
      internal::MacroscopicForceCorrection< LatticeModel_T >::apply( it, latticeModel, momentumDensity );
   }

   template< typename PdfField_T >
   static void get( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::getMomentumDensityD3Q19< LatticeModel_T, PdfField_T >( momentumDensity, pdf, x, y, z );
      internal::MacroscopicForceCorrection< LatticeModel_T >::apply( x, y, z, latticeModel, momentumDensity );
   }
};



} // namespace lbm
} // namespace walberla
