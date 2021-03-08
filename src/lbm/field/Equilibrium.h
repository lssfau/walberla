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
//! \file Equilibrium.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include "stencil/D3Q19.h"

#include <type_traits>


// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

//////////////////////////////////////////
// set equilibrium distribution (x,y,z) //
//////////////////////////////////////////

template< typename LatticeModel_T, class Enable = void >
struct Equilibrium
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::Equilibrium' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                             ! LatticeModel_T::compressible &&
                                                             LatticeModel_T::equilibriumAccuracyOrder == 2
                                                             >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                "Only works with incompressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t dir_independent = (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         it[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      real_t & xyz0 = pdf(x,y,z,0);
      const real_t dir_independent = (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         pdf.getF( &xyz0, d.toIdx() ) = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }
   }
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                             LatticeModel_T::equilibriumAccuracyOrder == 1
                                                             >::type >
{
   static_assert( LatticeModel_T::compressible == false,         "Only works with incompressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 1, "Only works for lattice models that require the equilibrium distribution to be order 1 accurate!" );

   // Right now there is now specialization for D3Q19 -> for all stencils this generic version is chosen

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t dir_independent = rho - real_t(1.0);
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         it[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel );
      }
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      real_t & xyz0 = pdf(x,y,z,0);
      const real_t dir_independent = rho - real_t(1.0);
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         pdf.getF( &xyz0, d.toIdx() ) = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel );
      }
   }
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                             LatticeModel_T::compressible &&
                                                             LatticeModel_T::equilibriumAccuracyOrder == 2
                                                             >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                         "Only works with compressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t dir_independent = real_t(1.0) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         it[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      real_t & xyz0 = pdf(x,y,z,0);
      const real_t dir_independent = real_t(1.0) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         pdf.getF( &xyz0, d.toIdx() ) = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }
   }
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
	                                                          LatticeModel_T::equilibriumAccuracyOrder == 1
	                                                          >::type >
{
   static_assert( LatticeModel_T::compressible,                  "Only works with compressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 1, "Only works for lattice models that require the equilibrium distribution to be order 1 accurate!" );

   // Right now there is now specialization for D3Q19 -> for all stencils this generic version is chosen

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         it[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( real_t(1.0) + real_t(3.0)*vel );
      }
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      real_t & xyz0 = pdf(x,y,z,0);
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         pdf.getF( &xyz0, d.toIdx() ) = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( real_t(1.0) + real_t(3.0)*vel );
      }
   }
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                             ! LatticeModel_T::compressible &&
                                                             LatticeModel_T::equilibriumAccuracyOrder == 2
                                                             >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                       "Only works with incompressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using Stencil = typename LatticeModel_T::Stencil;

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;

      const real_t velXX = velocity[0] * velocity[0];
      const real_t velYY = velocity[1] * velocity[1];
      const real_t velZZ = velocity[2] * velocity[2];

      const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) * (rho - real_t(1.0)) - real_t(0.5) * ( velXX + velYY + velZZ );

      it[ Stencil::idx[C] ] = dir_indep_trm;

      const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
      const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
      const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

      const real_t w1 = real_t(3.0) / real_t(18.0);

      it[ Stencil::idx[E] ] = w1 * ( vel_trm_E_W + velocity[0] );
      it[ Stencil::idx[W] ] = w1 * ( vel_trm_E_W - velocity[0] );
      it[ Stencil::idx[N] ] = w1 * ( vel_trm_N_S + velocity[1] );
      it[ Stencil::idx[S] ] = w1 * ( vel_trm_N_S - velocity[1] );
      it[ Stencil::idx[T] ] = w1 * ( vel_trm_T_B + velocity[2] );
      it[ Stencil::idx[B] ] = w1 * ( vel_trm_T_B - velocity[2] );

      const real_t velXmY = velocity[0] - velocity[1];
      const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

      const real_t w2 = real_t(3.0) / real_t(36.0);

      it[ Stencil::idx[NW] ] = w2 * ( vel_trm_NW_SE - velXmY );
      it[ Stencil::idx[SE] ] = w2 * ( vel_trm_NW_SE + velXmY );

      const real_t velXpY = velocity[0] + velocity[1];
      const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

      it[ Stencil::idx[NE] ] = w2 * ( vel_trm_NE_SW + velXpY );
      it[ Stencil::idx[SW] ] = w2 * ( vel_trm_NE_SW - velXpY );

      const real_t velXmZ = velocity[0] - velocity[2];
      const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

      it[ Stencil::idx[TW] ] = w2 * ( vel_trm_TW_BE - velXmZ );
      it[ Stencil::idx[BE] ] = w2 * ( vel_trm_TW_BE + velXmZ );

      const real_t velXpZ = velocity[0] + velocity[2];
      const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

      it[ Stencil::idx[TE] ] = w2 * ( vel_trm_TE_BW + velXpZ );
      it[ Stencil::idx[BW] ] = w2 * ( vel_trm_TE_BW - velXpZ );

      const real_t velYmZ = velocity[1] - velocity[2];
      const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

      it[ Stencil::idx[TS] ] = w2 * ( vel_trm_TS_BN - velYmZ );
      it[ Stencil::idx[BN] ] = w2 * ( vel_trm_TS_BN + velYmZ );

      const real_t velYpZ = velocity[1] + velocity[2];
      const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

      it[ Stencil::idx[TN] ] = w2 * ( vel_trm_TN_BS + velYpZ );
      it[ Stencil::idx[BS] ] = w2 * ( vel_trm_TN_BS - velYpZ );
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;

      real_t & xyz0 = pdf(x,y,z,0);

      const real_t velXX = velocity[0] * velocity[0];
      const real_t velYY = velocity[1] * velocity[1];
      const real_t velZZ = velocity[2] * velocity[2];

      const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) * (rho - real_t(1.0)) - real_t(0.5) * ( velXX + velYY + velZZ );

      pdf.getF( &xyz0, Stencil::idx[C] ) = dir_indep_trm;

      const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
      const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
      const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

      const real_t w1 = real_t(3.0) / real_t(18.0);

      pdf.getF( &xyz0, Stencil::idx[E] ) = w1 * ( vel_trm_E_W + velocity[0] );
      pdf.getF( &xyz0, Stencil::idx[W] ) = w1 * ( vel_trm_E_W - velocity[0] );
      pdf.getF( &xyz0, Stencil::idx[N] ) = w1 * ( vel_trm_N_S + velocity[1] );
      pdf.getF( &xyz0, Stencil::idx[S] ) = w1 * ( vel_trm_N_S - velocity[1] );
      pdf.getF( &xyz0, Stencil::idx[T] ) = w1 * ( vel_trm_T_B + velocity[2] );
      pdf.getF( &xyz0, Stencil::idx[B] ) = w1 * ( vel_trm_T_B - velocity[2] );

      const real_t velXmY = velocity[0] - velocity[1];
      const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

      const real_t w2 = real_t(3.0) / real_t(36.0);

      pdf.getF( &xyz0, Stencil::idx[NW] ) = w2 * ( vel_trm_NW_SE - velXmY );
      pdf.getF( &xyz0, Stencil::idx[SE] ) = w2 * ( vel_trm_NW_SE + velXmY );

      const real_t velXpY = velocity[0] + velocity[1];
      const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

      pdf.getF( &xyz0, Stencil::idx[NE] ) = w2 * ( vel_trm_NE_SW + velXpY );
      pdf.getF( &xyz0, Stencil::idx[SW] ) = w2 * ( vel_trm_NE_SW - velXpY );

      const real_t velXmZ = velocity[0] - velocity[2];
      const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

      pdf.getF( &xyz0, Stencil::idx[TW] ) = w2 * ( vel_trm_TW_BE - velXmZ );
      pdf.getF( &xyz0, Stencil::idx[BE] ) = w2 * ( vel_trm_TW_BE + velXmZ );

      const real_t velXpZ = velocity[0] + velocity[2];
      const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

      pdf.getF( &xyz0, Stencil::idx[TE] ) = w2 * ( vel_trm_TE_BW + velXpZ );
      pdf.getF( &xyz0, Stencil::idx[BW] ) = w2 * ( vel_trm_TE_BW - velXpZ );

      const real_t velYmZ = velocity[1] - velocity[2];
      const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

      pdf.getF( &xyz0, Stencil::idx[TS] ) = w2 * ( vel_trm_TS_BN - velYmZ );
      pdf.getF( &xyz0, Stencil::idx[BN] ) = w2 * ( vel_trm_TS_BN + velYmZ );

      const real_t velYpZ = velocity[1] + velocity[2];
      const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

      pdf.getF( &xyz0, Stencil::idx[TN] ) = w2 * ( vel_trm_TN_BS + velYpZ );
      pdf.getF( &xyz0, Stencil::idx[BS] ) = w2 * ( vel_trm_TN_BS - velYpZ );
   }
};



template< typename LatticeModel_T >
struct Equilibrium< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
	                                                          LatticeModel_T::compressible &&
	                                                          LatticeModel_T::equilibriumAccuracyOrder == 2
	                                                          >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                "Only works with compressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using Stencil = typename LatticeModel_T::Stencil;

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;

      const real_t velXX = velocity[0] * velocity[0];
      const real_t velYY = velocity[1] * velocity[1];
      const real_t velZZ = velocity[2] * velocity[2];

      const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) - real_t(0.5) * ( velXX + velYY + velZZ );

      it[ Stencil::idx[C] ] = rho * dir_indep_trm;

      const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
      const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
      const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

      const real_t w1_rho = rho * real_t(3.0) / real_t(18.0);

      it[ Stencil::idx[E] ] = w1_rho * ( vel_trm_E_W + velocity[0] );
      it[ Stencil::idx[W] ] = w1_rho * ( vel_trm_E_W - velocity[0] );
      it[ Stencil::idx[N] ] = w1_rho * ( vel_trm_N_S + velocity[1] );
      it[ Stencil::idx[S] ] = w1_rho * ( vel_trm_N_S - velocity[1] );
      it[ Stencil::idx[T] ] = w1_rho * ( vel_trm_T_B + velocity[2] );
      it[ Stencil::idx[B] ] = w1_rho * ( vel_trm_T_B - velocity[2] );

      const real_t velXmY = velocity[0] - velocity[1];
      const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

      const real_t w2_rho = rho * real_t(3.0) / real_t(36.0);

      it[ Stencil::idx[NW] ] = w2_rho * ( vel_trm_NW_SE - velXmY );
      it[ Stencil::idx[SE] ] = w2_rho * ( vel_trm_NW_SE + velXmY );

      const real_t velXpY = velocity[0] + velocity[1];
      const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

      it[ Stencil::idx[NE] ] = w2_rho * ( vel_trm_NE_SW + velXpY );
      it[ Stencil::idx[SW] ] = w2_rho * ( vel_trm_NE_SW - velXpY );

      const real_t velXmZ = velocity[0] - velocity[2];
      const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

      it[ Stencil::idx[TW] ] = w2_rho * ( vel_trm_TW_BE - velXmZ );
      it[ Stencil::idx[BE] ] = w2_rho * ( vel_trm_TW_BE + velXmZ );

      const real_t velXpZ = velocity[0] + velocity[2];
      const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

      it[ Stencil::idx[TE] ] = w2_rho * ( vel_trm_TE_BW + velXpZ );
      it[ Stencil::idx[BW] ] = w2_rho * ( vel_trm_TE_BW - velXpZ );

      const real_t velYmZ = velocity[1] - velocity[2];
      const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

      it[ Stencil::idx[TS] ] = w2_rho * ( vel_trm_TS_BN - velYmZ );
      it[ Stencil::idx[BN] ] = w2_rho * ( vel_trm_TS_BN + velYmZ );

      const real_t velYpZ = velocity[1] + velocity[2];
      const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

      it[ Stencil::idx[TN] ] = w2_rho * ( vel_trm_TN_BS + velYpZ );
      it[ Stencil::idx[BS] ] = w2_rho * ( vel_trm_TN_BS - velYpZ );
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;

      real_t & xyz0 = pdf(x,y,z,0);

      const real_t velXX = velocity[0] * velocity[0];
      const real_t velYY = velocity[1] * velocity[1];
      const real_t velZZ = velocity[2] * velocity[2];

      const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) - real_t(0.5) * ( velXX + velYY + velZZ );

      pdf.getF( &xyz0, Stencil::idx[C] ) = rho * dir_indep_trm;

      const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
      const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
      const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

      const real_t w1_rho = rho * real_t(3.0) / real_t(18.0);

      pdf.getF( &xyz0, Stencil::idx[E] ) = w1_rho * ( vel_trm_E_W + velocity[0] );
      pdf.getF( &xyz0, Stencil::idx[W] ) = w1_rho * ( vel_trm_E_W - velocity[0] );
      pdf.getF( &xyz0, Stencil::idx[N] ) = w1_rho * ( vel_trm_N_S + velocity[1] );
      pdf.getF( &xyz0, Stencil::idx[S] ) = w1_rho * ( vel_trm_N_S - velocity[1] );
      pdf.getF( &xyz0, Stencil::idx[T] ) = w1_rho * ( vel_trm_T_B + velocity[2] );
      pdf.getF( &xyz0, Stencil::idx[B] ) = w1_rho * ( vel_trm_T_B - velocity[2] );

      const real_t velXmY = velocity[0] - velocity[1];
      const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

      const real_t w2_rho = rho * real_t(3.0) / real_t(36.0);

      pdf.getF( &xyz0, Stencil::idx[NW] ) = w2_rho * ( vel_trm_NW_SE - velXmY );
      pdf.getF( &xyz0, Stencil::idx[SE] ) = w2_rho * ( vel_trm_NW_SE + velXmY );

      const real_t velXpY = velocity[0] + velocity[1];
      const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

      pdf.getF( &xyz0, Stencil::idx[NE] ) = w2_rho * ( vel_trm_NE_SW + velXpY );
      pdf.getF( &xyz0, Stencil::idx[SW] ) = w2_rho * ( vel_trm_NE_SW - velXpY );

      const real_t velXmZ = velocity[0] - velocity[2];
      const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

      pdf.getF( &xyz0, Stencil::idx[TW] ) = w2_rho * ( vel_trm_TW_BE - velXmZ );
      pdf.getF( &xyz0, Stencil::idx[BE] ) = w2_rho * ( vel_trm_TW_BE + velXmZ );

      const real_t velXpZ = velocity[0] + velocity[2];
      const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

      pdf.getF( &xyz0, Stencil::idx[TE] ) = w2_rho * ( vel_trm_TE_BW + velXpZ );
      pdf.getF( &xyz0, Stencil::idx[BW] ) = w2_rho * ( vel_trm_TE_BW - velXpZ );

      const real_t velYmZ = velocity[1] - velocity[2];
      const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

      pdf.getF( &xyz0, Stencil::idx[TS] ) = w2_rho * ( vel_trm_TS_BN - velYmZ );
      pdf.getF( &xyz0, Stencil::idx[BN] ) = w2_rho * ( vel_trm_TS_BN + velYmZ );

      const real_t velYpZ = velocity[1] + velocity[2];
      const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

      pdf.getF( &xyz0, Stencil::idx[TN] ) = w2_rho * ( vel_trm_TN_BS + velYpZ );
      pdf.getF( &xyz0, Stencil::idx[BS] ) = w2_rho * ( vel_trm_TN_BS - velYpZ );
   }
};



//////////////////////////////////////////
// set equilibrium distribution (range) //
//////////////////////////////////////////

template< typename LatticeModel_T, typename FieldIteratorXYZ, class Enable = void >
struct EquilibriumRange
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::EquilibriumRange' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct EquilibriumRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                                    LatticeModel_T::equilibriumAccuracyOrder == 2
                                                                                    >::type >
{
   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > value( stencil::NR_OF_DIRECTIONS );

      const real_t dir_independent = (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel =  real_c(d.cx()) * velocity[0] +  real_c(d.cy()) * velocity[1] +  real_c(d.cz()) * velocity[2];
         value[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }

      for( auto cell = begin; cell != end; ++cell )
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
            cell.getF( d.toIdx() ) = value[ d.toIdx() ];
   }
};

template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct EquilibriumRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                                    LatticeModel_T::equilibriumAccuracyOrder == 1
                                                                                    >::type >
{
   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > value( stencil::NR_OF_DIRECTIONS );

      const real_t dir_independent = rho - real_t(1.0);
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel =  real_c(d.cx()) * velocity[0] +  real_c(d.cy()) * velocity[1] +  real_c(d.cz()) * velocity[2];
         value[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * ( dir_independent + real_t(3.0)*vel );
      }

      for( auto cell = begin; cell != end; ++cell )
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
            cell.getF( d.toIdx() ) = value[ d.toIdx() ];
   }
};

template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct EquilibriumRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< LatticeModel_T::compressible &&
                                                                                    LatticeModel_T::equilibriumAccuracyOrder == 2
                                                                                    >::type >
{
   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > value( stencil::NR_OF_DIRECTIONS );

      const real_t dir_independent = real_t(1.0) - real_t(1.5) * velocity.sqrLength();
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         value[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( dir_independent + real_t(3.0)*vel + real_t(4.5)*vel*vel );
      }

      for( auto cell = begin; cell != end; ++cell )
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
            cell.getF( d.toIdx() ) = value[ d.toIdx() ];
   }
};

template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct EquilibriumRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< LatticeModel_T::compressible &&
                                                                                    LatticeModel_T::equilibriumAccuracyOrder == 1
                                                                                    >::type >
{
   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > value( stencil::NR_OF_DIRECTIONS );

      const real_t dir_independent = real_t(1.0);
      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
         value[ d.toIdx() ] = real_c(LatticeModel_T::w[ d.toIdx() ]) * rho * ( dir_independent + real_t(3.0)*vel );
      }

      for( auto cell = begin; cell != end; ++cell )
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
            cell.getF( d.toIdx() ) = value[ d.toIdx() ];
   }
};



} // namespace lbm
} // namespace walberla
