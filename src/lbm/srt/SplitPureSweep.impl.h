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
//! \file SplitPureSweep.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================



#ifndef __IBMCPP__

#pragma once

#include "lbm/lattice_model/LatticeModelBase.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/Streaming.h"
#include "lbm/sweeps/SweepBase.h"
#include "lbm/IntelCompilerOptimization.h"

#include "field/iterators/IteratorMacros.h"

#include <type_traits>


namespace walberla {
namespace lbm {

///////////////////////////////////////////////////////
// Available SRT implementations:                    //
//                                                   //
// There are no generic (D*Q*) versions!             //
//                                                   //
// Optimized D3Q19 implementation:                   //
//                     incompressible | compressible //
//          no forces:       x               x       //
///////////////////////////////////////////////////////


///////////////////////////////
// Specialization for D3Q19: //
// - incompressible          //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                               LatticeModel_T::CollisionModel::constant &&
                                                               std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                               ! LatticeModel_T::compressible &&
                                                               std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                  >::type > :
   public SweepBase< LatticeModel_T >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename SweepBase<LatticeModel_T>::PdfField_T  PdfField_T;
   typedef typename LatticeModel_T::Stencil                Stencil;

   // block has NO dst pdf field
   SplitPureSweep( const BlockDataID & pdfField ) :
      SweepBase<LatticeModel_T>( pdfField ) {}

   // every block has a dedicated dst pdf field
   SplitPureSweep( const BlockDataID & src, const BlockDataID & dst ) :
      SweepBase<LatticeModel_T>( src, dst ) {}

   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              ! LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
   >::operator()( IBlock * const block )
{
   PdfField_T * src( nullptr );
   PdfField_T * dst( nullptr );

   this->getFields( block, src, dst );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), 1 );

   // constants used during stream/collide

   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT dir_indep_trm = new real_t[ uint_c( xSize ) ];

   if( src->layout() == field::fzyx && dst->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         real_t * WALBERLA_RESTRICT pNE = &src->get(-1, y-1, z  , Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT pN  = &src->get(0 , y-1, z  , Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT pNW = &src->get(+1, y-1, z  , Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT pW  = &src->get(+1, y  , z  , Stencil::idx[W]);
         real_t * WALBERLA_RESTRICT pSW = &src->get(+1, y+1, z  , Stencil::idx[SW]);
         real_t * WALBERLA_RESTRICT pS  = &src->get(0 , y+1, z  , Stencil::idx[S]);
         real_t * WALBERLA_RESTRICT pSE = &src->get(-1, y+1, z  , Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT pE  = &src->get(-1, y  , z  , Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT pT  = &src->get(0 , y  , z-1, Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT pTE = &src->get(-1, y  , z-1, Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT pTN = &src->get(0 , y-1, z-1, Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT pTW = &src->get(+1, y  , z-1, Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT pTS = &src->get(0 , y+1, z-1, Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT pB  = &src->get(0 , y  , z+1, Stencil::idx[B]);
         real_t * WALBERLA_RESTRICT pBE = &src->get(-1, y  , z+1, Stencil::idx[BE]);
         real_t * WALBERLA_RESTRICT pBN = &src->get(0 , y-1, z+1, Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT pBW = &src->get(+1, y  , z+1, Stencil::idx[BW]);
         real_t * WALBERLA_RESTRICT pBS = &src->get(0 , y+1, z+1, Stencil::idx[BS]);
         real_t * WALBERLA_RESTRICT pC  = &src->get(0 , y  , z  , Stencil::idx[C]);

         real_t * WALBERLA_RESTRICT dC = &dst->get(0,y,z,Stencil::idx[C]);

         X_LOOP
         (
            const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
            const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
            const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

            const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;

            velX[x] = velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x];
            velY[x] = velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x];
            velZ[x] = velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x];

            dir_indep_trm[x] = one_third * rho - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            dC[x] = omega_trm * pC[x] + omega_w0 * dir_indep_trm[x];
         )

         real_t * WALBERLA_RESTRICT dNW = &dst->get(0,y,z,Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT dSE = &dst->get(0,y,z,Stencil::idx[SE]);

         X_LOOP
         (
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dNW[x] = omega_trm * pNW[x] + omega_w2 * ( vel_trm_NW_SE - vel );
            dSE[x] = omega_trm * pSE[x] + omega_w2 * ( vel_trm_NW_SE + vel );
         )

         real_t * WALBERLA_RESTRICT dNE = &dst->get(0,y,z,Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT dSW = &dst->get(0,y,z,Stencil::idx[SW]);

         X_LOOP
         (
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dNE[x] = omega_trm * pNE[x] + omega_w2 * ( vel_trm_NE_SW + vel );
            dSW[x] = omega_trm * pSW[x] + omega_w2 * ( vel_trm_NE_SW - vel );
         )

         real_t * WALBERLA_RESTRICT dTW = &dst->get(0,y,z,Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT dBE = &dst->get(0,y,z,Stencil::idx[BE]);

         X_LOOP
         (
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTW[x] = omega_trm * pTW[x] + omega_w2 * ( vel_trm_TW_BE - vel );
            dBE[x] = omega_trm * pBE[x] + omega_w2 * ( vel_trm_TW_BE + vel );
         )

         real_t * WALBERLA_RESTRICT dTE = &dst->get(0,y,z,Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT dBW = &dst->get(0,y,z,Stencil::idx[BW]);

         X_LOOP
         (
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTE[x] = omega_trm * pTE[x] + omega_w2 * ( vel_trm_TE_BW + vel );
            dBW[x] = omega_trm * pBW[x] + omega_w2 * ( vel_trm_TE_BW - vel );
         )

         real_t * WALBERLA_RESTRICT dTS = &dst->get(0,y,z,Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT dBN = &dst->get(0,y,z,Stencil::idx[BN]);

         X_LOOP
         (
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTS[x] = omega_trm * pTS[x] + omega_w2 * ( vel_trm_TS_BN - vel );
            dBN[x] = omega_trm * pBN[x] + omega_w2 * ( vel_trm_TS_BN + vel );
         )

         real_t * WALBERLA_RESTRICT dTN = &dst->get(0,y,z,Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT dBS = &dst->get(0,y,z,Stencil::idx[BS]);

         X_LOOP
         (
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTN[x] = omega_trm * pTN[x] + omega_w2 * ( vel_trm_TN_BS + vel );
            dBS[x] = omega_trm * pBS[x] + omega_w2 * ( vel_trm_TN_BS - vel );
         )

         real_t * WALBERLA_RESTRICT dN = &dst->get(0,y,z,Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT dS = &dst->get(0,y,z,Stencil::idx[S]);

         X_LOOP
         (
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            dN[x] = omega_trm * pN[x] + omega_w1 * ( vel_trm_N_S + velY[x] );
            dS[x] = omega_trm * pS[x] + omega_w1 * ( vel_trm_N_S - velY[x] );
         )

         real_t * WALBERLA_RESTRICT dE = &dst->get(0,y,z,Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT dW = &dst->get(0,y,z,Stencil::idx[W]);

         X_LOOP
         (
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            dE[x] = omega_trm * pE[x] + omega_w1 * ( vel_trm_E_W + velX[x] );
            dW[x] = omega_trm * pW[x] + omega_w1 * ( vel_trm_E_W - velX[x] );
         )

         real_t * WALBERLA_RESTRICT dT = &dst->get(0,y,z,Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT dB = &dst->get(0,y,z,Stencil::idx[B]);

         X_LOOP
         (
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            dT[x] = omega_trm * pT[x] + omega_w1 * ( vel_trm_T_B + velZ[x] );
            dB[x] = omega_trm * pB[x] + omega_w1 * ( vel_trm_T_B - velZ[x] );
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t dd_tmp_NE = src->get(x-1, y-1, z  , Stencil::idx[NE]);
            const real_t dd_tmp_N  = src->get(x  , y-1, z  , Stencil::idx[N]);
            const real_t dd_tmp_NW = src->get(x+1, y-1, z  , Stencil::idx[NW]);
            const real_t dd_tmp_W  = src->get(x+1, y  , z  , Stencil::idx[W]);
            const real_t dd_tmp_SW = src->get(x+1, y+1, z  , Stencil::idx[SW]);
            const real_t dd_tmp_S  = src->get(x  , y+1, z  , Stencil::idx[S]);
            const real_t dd_tmp_SE = src->get(x-1, y+1, z  , Stencil::idx[SE]);
            const real_t dd_tmp_E  = src->get(x-1, y  , z  , Stencil::idx[E]);
            const real_t dd_tmp_T  = src->get(x  , y  , z-1, Stencil::idx[T]);
            const real_t dd_tmp_TE = src->get(x-1, y  , z-1, Stencil::idx[TE]);
            const real_t dd_tmp_TN = src->get(x  , y-1, z-1, Stencil::idx[TN]);
            const real_t dd_tmp_TW = src->get(x+1, y  , z-1, Stencil::idx[TW]);
            const real_t dd_tmp_TS = src->get(x  , y+1, z-1, Stencil::idx[TS]);
            const real_t dd_tmp_B  = src->get(x  , y  , z+1, Stencil::idx[B]);
            const real_t dd_tmp_BE = src->get(x-1, y  , z+1, Stencil::idx[BE]);
            const real_t dd_tmp_BN = src->get(x  , y-1, z+1, Stencil::idx[BN]);
            const real_t dd_tmp_BW = src->get(x+1, y  , z+1, Stencil::idx[BW]);
            const real_t dd_tmp_BS = src->get(x  , y+1, z+1, Stencil::idx[BS]);
            const real_t dd_tmp_C  = src->get(x  , y  , z  , Stencil::idx[C]);

            const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
            const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
            const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

            const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

            velX[x] = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
            velY[x] = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
            velZ[x] = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

            dir_indep_trm[x] = one_third * rho - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            dst->get(x,y,z,Stencil::idx[C]) = omega_trm * dd_tmp_C + omega_w0 * dir_indep_trm[x];
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[NW]) = omega_trm * src->get(x+1, y-1, z, Stencil::idx[NW]) + omega_w2 * ( vel_trm_NW_SE - vel );
            dst->get(x,y,z,Stencil::idx[SE]) = omega_trm * src->get(x-1, y+1, z, Stencil::idx[SE]) + omega_w2 * ( vel_trm_NW_SE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[NE]) = omega_trm * src->get(x-1, y-1, z, Stencil::idx[NE]) + omega_w2 * ( vel_trm_NE_SW + vel );
            dst->get(x,y,z,Stencil::idx[SW]) = omega_trm * src->get(x+1, y+1, z, Stencil::idx[SW]) + omega_w2 * ( vel_trm_NE_SW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TW]) = omega_trm * src->get(x+1, y, z-1, Stencil::idx[TW]) + omega_w2 * ( vel_trm_TW_BE - vel );
            dst->get(x,y,z,Stencil::idx[BE]) = omega_trm * src->get(x-1, y, z+1, Stencil::idx[BE]) + omega_w2 * ( vel_trm_TW_BE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TE]) = omega_trm * src->get(x-1, y, z-1, Stencil::idx[TE]) + omega_w2 * ( vel_trm_TE_BW + vel );
            dst->get(x,y,z,Stencil::idx[BW]) = omega_trm * src->get(x+1, y, z+1, Stencil::idx[BW]) + omega_w2 * ( vel_trm_TE_BW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TS]) = omega_trm * src->get(x, y+1, z-1, Stencil::idx[TS]) + omega_w2 * ( vel_trm_TS_BN - vel );
            dst->get(x,y,z,Stencil::idx[BN]) = omega_trm * src->get(x, y-1, z+1, Stencil::idx[BN]) + omega_w2 * ( vel_trm_TS_BN + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TN]) = omega_trm * src->get(x, y-1, z-1, Stencil::idx[TN]) + omega_w2 * ( vel_trm_TN_BS + vel );
            dst->get(x,y,z,Stencil::idx[BS]) = omega_trm * src->get(x, y+1, z+1, Stencil::idx[BS]) + omega_w2 * ( vel_trm_TN_BS - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            dst->get(x,y,z,Stencil::idx[N]) = omega_trm * src->get(x, y-1, z, Stencil::idx[N]) + omega_w1 * ( vel_trm_N_S + velY[x] );
            dst->get(x,y,z,Stencil::idx[S]) = omega_trm * src->get(x, y+1, z, Stencil::idx[S]) + omega_w1 * ( vel_trm_N_S - velY[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            dst->get(x,y,z,Stencil::idx[E]) = omega_trm * src->get(x-1, y, z, Stencil::idx[E]) + omega_w1 * ( vel_trm_E_W + velX[x] );
            dst->get(x,y,z,Stencil::idx[W]) = omega_trm * src->get(x+1, y, z, Stencil::idx[W]) + omega_w1 * ( vel_trm_E_W - velX[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            dst->get(x,y,z,Stencil::idx[T]) = omega_trm * src->get(x, y, z-1, Stencil::idx[T]) + omega_w1 * ( vel_trm_T_B + velZ[x] );
            dst->get(x,y,z,Stencil::idx[B]) = omega_trm * src->get(x, y, z+1, Stencil::idx[B]) + omega_w1 * ( vel_trm_T_B - velZ[x] );
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] dir_indep_trm;

#ifdef _OPENMP
   }
#endif

   src->swapDataPointers( dst );
}



template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              ! LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
   >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( nullptr );
   PdfField_T * dst( nullptr );

   this->getFields( block, src, dst );

   StreamEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}



template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              ! LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
#ifdef NDEBUG
   >::collide( IBlock * const block, const uint_t /*numberOfGhostLayersToInclude*/ )
#else
   >::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
#endif
{
   WALBERLA_ASSERT_EQUAL( numberOfGhostLayersToInclude, uint_t(0) ); // the implementation right now doesn't support inclusion of ghost layers in collide step!

   PdfField_T * src = this->getSrcField( block );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   // constants used during stream/collide

   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT dir_indep_trm = new real_t[ uint_c( xSize ) ];

   if( src->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         real_t * WALBERLA_RESTRICT pC  = &src->get( 0, y, z, Stencil::idx[C]);
         real_t * WALBERLA_RESTRICT pN  = &src->get( 0, y, z, Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT pS  = &src->get( 0, y, z, Stencil::idx[S]);
         real_t * WALBERLA_RESTRICT pW  = &src->get( 0, y, z, Stencil::idx[W]);
         real_t * WALBERLA_RESTRICT pE  = &src->get( 0, y, z, Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT pT  = &src->get( 0, y, z, Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT pB  = &src->get( 0, y, z, Stencil::idx[B]);
         real_t * WALBERLA_RESTRICT pNW = &src->get( 0, y, z, Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT pNE = &src->get( 0, y, z, Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT pSW = &src->get( 0, y, z, Stencil::idx[SW]);
         real_t * WALBERLA_RESTRICT pSE = &src->get( 0, y, z, Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT pTN = &src->get( 0, y, z, Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT pTS = &src->get( 0, y, z, Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT pTW = &src->get( 0, y, z, Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT pTE = &src->get( 0, y, z, Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT pBN = &src->get( 0, y, z, Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT pBS = &src->get( 0, y, z, Stencil::idx[BS]);
         real_t * WALBERLA_RESTRICT pBW = &src->get( 0, y, z, Stencil::idx[BW]);
         real_t * WALBERLA_RESTRICT pBE = &src->get( 0, y, z, Stencil::idx[BE]);

         X_LOOP
         (
            const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
            const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
            const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

            const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;

            velX[x] = velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x];
            velY[x] = velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x];
            velZ[x] = velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x];

            dir_indep_trm[x] = one_third * rho - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            pC[x] = omega_trm * pC[x] + omega_w0 * dir_indep_trm[x];
         )

         X_LOOP
         (
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pNW[x] = omega_trm * pNW[x] + omega_w2 * ( vel_trm_NW_SE - vel );
            pSE[x] = omega_trm * pSE[x] + omega_w2 * ( vel_trm_NW_SE + vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pNE[x] = omega_trm * pNE[x] + omega_w2 * ( vel_trm_NE_SW + vel );
            pSW[x] = omega_trm * pSW[x] + omega_w2 * ( vel_trm_NE_SW - vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTW[x] = omega_trm * pTW[x] + omega_w2 * ( vel_trm_TW_BE - vel );
            pBE[x] = omega_trm * pBE[x] + omega_w2 * ( vel_trm_TW_BE + vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTE[x] = omega_trm * pTE[x] + omega_w2 * ( vel_trm_TE_BW + vel );
            pBW[x] = omega_trm * pBW[x] + omega_w2 * ( vel_trm_TE_BW - vel );
         )

         X_LOOP
         (
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTS[x] = omega_trm * pTS[x] + omega_w2 * ( vel_trm_TS_BN - vel );
            pBN[x] = omega_trm * pBN[x] + omega_w2 * ( vel_trm_TS_BN + vel );
         )

         X_LOOP
         (
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTN[x] = omega_trm * pTN[x] + omega_w2 * ( vel_trm_TN_BS + vel );
            pBS[x] = omega_trm * pBS[x] + omega_w2 * ( vel_trm_TN_BS - vel );
         )

         X_LOOP
         (
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            pN[x] = omega_trm * pN[x] + omega_w1 * ( vel_trm_N_S + velY[x] );
            pS[x] = omega_trm * pS[x] + omega_w1 * ( vel_trm_N_S - velY[x] );
         )

         X_LOOP
         (
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            pE[x] = omega_trm * pE[x] + omega_w1 * ( vel_trm_E_W + velX[x] );
            pW[x] = omega_trm * pW[x] + omega_w1 * ( vel_trm_E_W - velX[x] );
         )

         X_LOOP
         (
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            pT[x] = omega_trm * pT[x] + omega_w1 * ( vel_trm_T_B + velZ[x] );
            pB[x] = omega_trm * pB[x] + omega_w1 * ( vel_trm_T_B - velZ[x] );
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t dd_tmp_C  = src->get( x, y, z, Stencil::idx[C]  );
            const real_t dd_tmp_N  = src->get( x, y, z, Stencil::idx[N]  );
            const real_t dd_tmp_S  = src->get( x, y, z, Stencil::idx[S]  );
            const real_t dd_tmp_W  = src->get( x, y, z, Stencil::idx[W]  );
            const real_t dd_tmp_E  = src->get( x, y, z, Stencil::idx[E]  );
            const real_t dd_tmp_T  = src->get( x, y, z, Stencil::idx[T]  );
            const real_t dd_tmp_B  = src->get( x, y, z, Stencil::idx[B]  );
            const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW] );
            const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE] );
            const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW] );
            const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE] );
            const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN] );
            const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS] );
            const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW] );
            const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE] );
            const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN] );
            const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS] );
            const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW] );
            const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE] );

            const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
            const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
            const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

            const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

            velX[x] = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
            velY[x] = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
            velZ[x] = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

            dir_indep_trm[x] = one_third * rho - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            src->get(x,y,z,Stencil::idx[C]) = omega_trm * dd_tmp_C + omega_w0 * dir_indep_trm[x];
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[NW]) = omega_trm * src->get(x,y,z,Stencil::idx[NW]) + omega_w2 * ( vel_trm_NW_SE - vel );
            src->get(x,y,z,Stencil::idx[SE]) = omega_trm * src->get(x,y,z,Stencil::idx[SE]) + omega_w2 * ( vel_trm_NW_SE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[NE]) = omega_trm * src->get(x,y,z,Stencil::idx[NE]) + omega_w2 * ( vel_trm_NE_SW + vel );
            src->get(x,y,z,Stencil::idx[SW]) = omega_trm * src->get(x,y,z,Stencil::idx[SW]) + omega_w2 * ( vel_trm_NE_SW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TW]) = omega_trm * src->get(x,y,z,Stencil::idx[TW]) + omega_w2 * ( vel_trm_TW_BE - vel );
            src->get(x,y,z,Stencil::idx[BE]) = omega_trm * src->get(x,y,z,Stencil::idx[BE]) + omega_w2 * ( vel_trm_TW_BE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TE]) = omega_trm * src->get(x,y,z,Stencil::idx[TE]) + omega_w2 * ( vel_trm_TE_BW + vel );
            src->get(x,y,z,Stencil::idx[BW]) = omega_trm * src->get(x,y,z,Stencil::idx[BW]) + omega_w2 * ( vel_trm_TE_BW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TS]) = omega_trm * src->get(x,y,z,Stencil::idx[TS]) + omega_w2 * ( vel_trm_TS_BN - vel );
            src->get(x,y,z,Stencil::idx[BN]) = omega_trm * src->get(x,y,z,Stencil::idx[BN]) + omega_w2 * ( vel_trm_TS_BN + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TN]) = omega_trm * src->get(x,y,z,Stencil::idx[TN]) + omega_w2 * ( vel_trm_TN_BS + vel );
            src->get(x,y,z,Stencil::idx[BS]) = omega_trm * src->get(x,y,z,Stencil::idx[BS]) + omega_w2 * ( vel_trm_TN_BS - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            src->get(x,y,z,Stencil::idx[N]) = omega_trm * src->get(x,y,z,Stencil::idx[N]) + omega_w1 * ( vel_trm_N_S + velY[x] );
            src->get(x,y,z,Stencil::idx[S]) = omega_trm * src->get(x,y,z,Stencil::idx[S]) + omega_w1 * ( vel_trm_N_S - velY[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            src->get(x,y,z,Stencil::idx[E]) = omega_trm * src->get(x,y,z,Stencil::idx[E]) + omega_w1 * ( vel_trm_E_W + velX[x] );
            src->get(x,y,z,Stencil::idx[W]) = omega_trm * src->get(x,y,z,Stencil::idx[W]) + omega_w1 * ( vel_trm_E_W - velX[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            src->get(x,y,z,Stencil::idx[T]) = omega_trm * src->get(x,y,z,Stencil::idx[T]) + omega_w1 * ( vel_trm_T_B + velZ[x] );
            src->get(x,y,z,Stencil::idx[B]) = omega_trm * src->get(x,y,z,Stencil::idx[B]) + omega_w1 * ( vel_trm_T_B - velZ[x] );
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] dir_indep_trm;

#ifdef _OPENMP
   }
#endif
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible          //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                               LatticeModel_T::CollisionModel::constant &&
                                                               std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                               LatticeModel_T::compressible &&
                                                               std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                  >::type > :
   public SweepBase< LatticeModel_T >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename SweepBase<LatticeModel_T>::PdfField_T  PdfField_T;
   typedef typename LatticeModel_T::Stencil                Stencil;

   // block has NO dst pdf field
   SplitPureSweep( const BlockDataID & pdfField ) :
      SweepBase<LatticeModel_T>( pdfField ) {}

   // every block has a dedicated dst pdf field
   SplitPureSweep( const BlockDataID & src, const BlockDataID & dst ) :
      SweepBase<LatticeModel_T>( src, dst ) {}

   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
   >::operator()( IBlock * const block )
{
   PdfField_T * src( nullptr );
   PdfField_T * dst( nullptr );

   this->getFields( block, src, dst );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), 1 );

   // constants used during stream/collide

   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT rho = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT dir_indep_trm = new real_t[ uint_c( xSize ) ];

   if( src->layout() == field::fzyx && dst->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         real_t * WALBERLA_RESTRICT pNE = &src->get(-1, y-1, z  , Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT pN  = &src->get(0 , y-1, z  , Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT pNW = &src->get(+1, y-1, z  , Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT pW  = &src->get(+1, y  , z  , Stencil::idx[W]);
         real_t * WALBERLA_RESTRICT pSW = &src->get(+1, y+1, z  , Stencil::idx[SW]);
         real_t * WALBERLA_RESTRICT pS  = &src->get(0 , y+1, z  , Stencil::idx[S]);
         real_t * WALBERLA_RESTRICT pSE = &src->get(-1, y+1, z  , Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT pE  = &src->get(-1, y  , z  , Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT pT  = &src->get(0 , y  , z-1, Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT pTE = &src->get(-1, y  , z-1, Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT pTN = &src->get(0 , y-1, z-1, Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT pTW = &src->get(+1, y  , z-1, Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT pTS = &src->get(0 , y+1, z-1, Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT pB  = &src->get(0 , y  , z+1, Stencil::idx[B]);
         real_t * WALBERLA_RESTRICT pBE = &src->get(-1, y  , z+1, Stencil::idx[BE]);
         real_t * WALBERLA_RESTRICT pBN = &src->get(0 , y-1, z+1, Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT pBW = &src->get(+1, y  , z+1, Stencil::idx[BW]);
         real_t * WALBERLA_RESTRICT pBS = &src->get(0 , y+1, z+1, Stencil::idx[BS]);
         real_t * WALBERLA_RESTRICT pC  = &src->get(0 , y  , z  , Stencil::idx[C]);

         real_t * WALBERLA_RESTRICT dC = &dst->get(0,y,z,Stencil::idx[C]);

         X_LOOP
         (
              const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
              const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
              const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

              rho[x] = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;
              const real_t rho_inv = real_t(1) / rho[x];

              velX[x] = rho_inv * ( velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x] );
              velY[x] = rho_inv * ( velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x] );
              velZ[x] = rho_inv * ( velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x] );

              dir_indep_trm[x] = one_third - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

              dC[x] = omega_trm * pC[x] + omega_w0 * rho[x] * dir_indep_trm[x];
         )

         real_t * WALBERLA_RESTRICT dNW = &dst->get(0,y,z,Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT dSE = &dst->get(0,y,z,Stencil::idx[SE]);

         X_LOOP
         (
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dNW[x] = omega_trm * pNW[x] + omega_w2 * rho[x] * ( vel_trm_NW_SE - vel );
            dSE[x] = omega_trm * pSE[x] + omega_w2 * rho[x] * ( vel_trm_NW_SE + vel );
         )

         real_t * WALBERLA_RESTRICT dNE = &dst->get(0,y,z,Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT dSW = &dst->get(0,y,z,Stencil::idx[SW]);

         X_LOOP
         (
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dNE[x] = omega_trm * pNE[x] + omega_w2 * rho[x] * ( vel_trm_NE_SW + vel );
            dSW[x] = omega_trm * pSW[x] + omega_w2 * rho[x] * ( vel_trm_NE_SW - vel );
         )

         real_t * WALBERLA_RESTRICT dTW = &dst->get(0,y,z,Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT dBE = &dst->get(0,y,z,Stencil::idx[BE]);

         X_LOOP
         (
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTW[x] = omega_trm * pTW[x] + omega_w2 * rho[x] * ( vel_trm_TW_BE - vel );
            dBE[x] = omega_trm * pBE[x] + omega_w2 * rho[x] * ( vel_trm_TW_BE + vel );
         )

         real_t * WALBERLA_RESTRICT dTE = &dst->get(0,y,z,Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT dBW = &dst->get(0,y,z,Stencil::idx[BW]);


         X_LOOP
         (
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTE[x] = omega_trm * pTE[x] + omega_w2 * rho[x] * ( vel_trm_TE_BW + vel );
            dBW[x] = omega_trm * pBW[x] + omega_w2 * rho[x] * ( vel_trm_TE_BW - vel );
         )

         real_t * WALBERLA_RESTRICT dTS = &dst->get(0,y,z,Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT dBN = &dst->get(0,y,z,Stencil::idx[BN]);

         X_LOOP
         (
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTS[x] = omega_trm * pTS[x] + omega_w2 * rho[x] * ( vel_trm_TS_BN - vel );
            dBN[x] = omega_trm * pBN[x] + omega_w2 * rho[x] * ( vel_trm_TS_BN + vel );
         )

         real_t * WALBERLA_RESTRICT dTN = &dst->get(0,y,z,Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT dBS = &dst->get(0,y,z,Stencil::idx[BS]);

         X_LOOP
         (
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dTN[x] = omega_trm * pTN[x] + omega_w2 * rho[x] * ( vel_trm_TN_BS + vel );
            dBS[x] = omega_trm * pBS[x] + omega_w2 * rho[x] * ( vel_trm_TN_BS - vel );
         )

         real_t * WALBERLA_RESTRICT dN = &dst->get(0,y,z,Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT dS = &dst->get(0,y,z,Stencil::idx[S]);

         X_LOOP
         (
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            dN[x] = omega_trm * pN[x] + omega_w1 * rho[x] * ( vel_trm_N_S + velY[x] );
            dS[x] = omega_trm * pS[x] + omega_w1 * rho[x] * ( vel_trm_N_S - velY[x] );
         )

         real_t * WALBERLA_RESTRICT dE = &dst->get(0,y,z,Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT dW = &dst->get(0,y,z,Stencil::idx[W]);

         X_LOOP
         (
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            dE[x] = omega_trm * pE[x] + omega_w1 * rho[x] * ( vel_trm_E_W + velX[x] );
            dW[x] = omega_trm * pW[x] + omega_w1 * rho[x] * ( vel_trm_E_W - velX[x] );
         )

         real_t * WALBERLA_RESTRICT dT = &dst->get(0,y,z,Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT dB = &dst->get(0,y,z,Stencil::idx[B]);

         X_LOOP
         (
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            dT[x] = omega_trm * pT[x] + omega_w1 * rho[x] * ( vel_trm_T_B + velZ[x] );
            dB[x] = omega_trm * pB[x] + omega_w1 * rho[x] * ( vel_trm_T_B - velZ[x] );
         )

       ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t dd_tmp_NE = src->get(x-1, y-1, z  , Stencil::idx[NE]);
            const real_t dd_tmp_N  = src->get(x  , y-1, z  , Stencil::idx[N]);
            const real_t dd_tmp_NW = src->get(x+1, y-1, z  , Stencil::idx[NW]);
            const real_t dd_tmp_W  = src->get(x+1, y  , z  , Stencil::idx[W]);
            const real_t dd_tmp_SW = src->get(x+1, y+1, z  , Stencil::idx[SW]);
            const real_t dd_tmp_S  = src->get(x  , y+1, z  , Stencil::idx[S]);
            const real_t dd_tmp_SE = src->get(x-1, y+1, z  , Stencil::idx[SE]);
            const real_t dd_tmp_E  = src->get(x-1, y  , z  , Stencil::idx[E]);
            const real_t dd_tmp_T  = src->get(x  , y  , z-1, Stencil::idx[T]);
            const real_t dd_tmp_TE = src->get(x-1, y  , z-1, Stencil::idx[TE]);
            const real_t dd_tmp_TN = src->get(x  , y-1, z-1, Stencil::idx[TN]);
            const real_t dd_tmp_TW = src->get(x+1, y  , z-1, Stencil::idx[TW]);
            const real_t dd_tmp_TS = src->get(x  , y+1, z-1, Stencil::idx[TS]);
            const real_t dd_tmp_B  = src->get(x  , y  , z+1, Stencil::idx[B]);
            const real_t dd_tmp_BE = src->get(x-1, y  , z+1, Stencil::idx[BE]);
            const real_t dd_tmp_BN = src->get(x  , y-1, z+1, Stencil::idx[BN]);
            const real_t dd_tmp_BW = src->get(x+1, y  , z+1, Stencil::idx[BW]);
            const real_t dd_tmp_BS = src->get(x  , y+1, z+1, Stencil::idx[BS]);
            const real_t dd_tmp_C  = src->get(x  , y  , z  , Stencil::idx[C]);

            const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
            const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
            const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

            rho[x] = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;
            const real_t rho_inv = real_t(1) / rho[x];

            velX[x] = rho_inv * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
            velY[x] = rho_inv * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
            velZ[x] = rho_inv * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

            dir_indep_trm[x] = one_third - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            dst->get(x,y,z,Stencil::idx[C]) = omega_trm * dd_tmp_C + omega_w0 * rho[x] * dir_indep_trm[x];
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[NW]) = omega_trm * src->get(x+1, y-1, z, Stencil::idx[NW]) + omega_w2 * rho[x] * ( vel_trm_NW_SE - vel );
            dst->get(x,y,z,Stencil::idx[SE]) = omega_trm * src->get(x-1, y+1, z, Stencil::idx[SE]) + omega_w2 * rho[x] * ( vel_trm_NW_SE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[NE]) = omega_trm * src->get(x-1, y-1, z, Stencil::idx[NE]) + omega_w2 * rho[x] * ( vel_trm_NE_SW + vel );
            dst->get(x,y,z,Stencil::idx[SW]) = omega_trm * src->get(x+1, y+1, z, Stencil::idx[SW]) + omega_w2 * rho[x] * ( vel_trm_NE_SW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TW]) = omega_trm * src->get(x+1, y, z-1, Stencil::idx[TW]) + omega_w2 * rho[x] * ( vel_trm_TW_BE - vel );
            dst->get(x,y,z,Stencil::idx[BE]) = omega_trm * src->get(x-1, y, z+1, Stencil::idx[BE]) + omega_w2 * rho[x] * ( vel_trm_TW_BE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TE]) = omega_trm * src->get(x-1, y, z-1, Stencil::idx[TE]) + omega_w2 * rho[x] * ( vel_trm_TE_BW + vel );
            dst->get(x,y,z,Stencil::idx[BW]) = omega_trm * src->get(x+1, y, z+1, Stencil::idx[BW]) + omega_w2 * rho[x] * ( vel_trm_TE_BW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TS]) = omega_trm * src->get(x, y+1, z-1, Stencil::idx[TS]) + omega_w2 * rho[x] * ( vel_trm_TS_BN - vel );
            dst->get(x,y,z,Stencil::idx[BN]) = omega_trm * src->get(x, y-1, z+1, Stencil::idx[BN]) + omega_w2 * rho[x] * ( vel_trm_TS_BN + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            dst->get(x,y,z,Stencil::idx[TN]) = omega_trm * src->get(x, y-1, z-1, Stencil::idx[TN]) + omega_w2 * rho[x] * ( vel_trm_TN_BS + vel );
            dst->get(x,y,z,Stencil::idx[BS]) = omega_trm * src->get(x, y+1, z+1, Stencil::idx[BS]) + omega_w2 * rho[x] * ( vel_trm_TN_BS - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            dst->get(x,y,z,Stencil::idx[N]) = omega_trm * src->get(x, y-1, z, Stencil::idx[N]) + omega_w1 * rho[x] * ( vel_trm_N_S + velY[x] );
            dst->get(x,y,z,Stencil::idx[S]) = omega_trm * src->get(x, y+1, z, Stencil::idx[S]) + omega_w1 * rho[x] * ( vel_trm_N_S - velY[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            dst->get(x,y,z,Stencil::idx[E]) = omega_trm * src->get(x-1, y, z, Stencil::idx[E]) + omega_w1 * rho[x] * ( vel_trm_E_W + velX[x] );
            dst->get(x,y,z,Stencil::idx[W]) = omega_trm * src->get(x+1, y, z, Stencil::idx[W]) + omega_w1 * rho[x] * ( vel_trm_E_W - velX[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            dst->get(x,y,z,Stencil::idx[T]) = omega_trm * src->get(x, y, z-1, Stencil::idx[T]) + omega_w1 * rho[x] * ( vel_trm_T_B + velZ[x] );
            dst->get(x,y,z,Stencil::idx[B]) = omega_trm * src->get(x, y, z+1, Stencil::idx[B]) + omega_w1 * rho[x] * ( vel_trm_T_B - velZ[x] );
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] rho;
   delete[] dir_indep_trm;

#ifdef _OPENMP
   }
#endif

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
   >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( nullptr );
   PdfField_T * dst( nullptr );

   this->getFields( block, src, dst );

   StreamEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T >
void SplitPureSweep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                              LatticeModel_T::CollisionModel::constant &&
                                                              std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                              LatticeModel_T::compressible &&
                                                              std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                 >::type
#ifdef NDEBUG                                                                    
   >::collide( IBlock * const block, const uint_t /*numberOfGhostLayersToInclude*/ )
#else
   >::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
#endif
{
   WALBERLA_ASSERT_EQUAL( numberOfGhostLayersToInclude, uint_t(0) ); // the implementation right now doesn't support inclusion of ghost layers in collide step!

   PdfField_T * src = this->getSrcField( block );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   // constants used during stream/collide

   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c( xSize ) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT rho = new real_t[ uint_c( xSize ) ];

   real_t * WALBERLA_RESTRICT dir_indep_trm = new real_t[ uint_c( xSize ) ];

   if( src->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         real_t * WALBERLA_RESTRICT pC  = &src->get( 0, y, z, Stencil::idx[C]);
         real_t * WALBERLA_RESTRICT pN  = &src->get( 0, y, z, Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT pS  = &src->get( 0, y, z, Stencil::idx[S]);
         real_t * WALBERLA_RESTRICT pW  = &src->get( 0, y, z, Stencil::idx[W]);
         real_t * WALBERLA_RESTRICT pE  = &src->get( 0, y, z, Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT pT  = &src->get( 0, y, z, Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT pB  = &src->get( 0, y, z, Stencil::idx[B]);
         real_t * WALBERLA_RESTRICT pNW = &src->get( 0, y, z, Stencil::idx[NW]);
         real_t * WALBERLA_RESTRICT pNE = &src->get( 0, y, z, Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT pSW = &src->get( 0, y, z, Stencil::idx[SW]);
         real_t * WALBERLA_RESTRICT pSE = &src->get( 0, y, z, Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT pTN = &src->get( 0, y, z, Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT pTS = &src->get( 0, y, z, Stencil::idx[TS]);
         real_t * WALBERLA_RESTRICT pTW = &src->get( 0, y, z, Stencil::idx[TW]);
         real_t * WALBERLA_RESTRICT pTE = &src->get( 0, y, z, Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT pBN = &src->get( 0, y, z, Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT pBS = &src->get( 0, y, z, Stencil::idx[BS]);
         real_t * WALBERLA_RESTRICT pBW = &src->get( 0, y, z, Stencil::idx[BW]);
         real_t * WALBERLA_RESTRICT pBE = &src->get( 0, y, z, Stencil::idx[BE]);

         X_LOOP
         (
            const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
            const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
            const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

            rho[x] = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;
            const real_t rho_inv = real_t(1) / rho[x];

            velX[x] = rho_inv * ( velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x] );
            velY[x] = rho_inv * ( velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x] );
            velZ[x] = rho_inv * ( velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x] );

            dir_indep_trm[x] = one_third - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            pC[x] = omega_trm * pC[x] + omega_w0 * rho[x] * dir_indep_trm[x];
         )

         X_LOOP
         (
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pNW[x] = omega_trm * pNW[x] + omega_w2 * rho[x] * ( vel_trm_NW_SE - vel );
            pSE[x] = omega_trm * pSE[x] + omega_w2 * rho[x] * ( vel_trm_NW_SE + vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pNE[x] = omega_trm * pNE[x] + omega_w2 * rho[x] * ( vel_trm_NE_SW + vel );
            pSW[x] = omega_trm * pSW[x] + omega_w2 * rho[x] * ( vel_trm_NE_SW - vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTW[x] = omega_trm * pTW[x] + omega_w2 * rho[x] * ( vel_trm_TW_BE - vel );
            pBE[x] = omega_trm * pBE[x] + omega_w2 * rho[x] * ( vel_trm_TW_BE + vel );
         )

         X_LOOP
         (
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTE[x] = omega_trm * pTE[x] + omega_w2 * rho[x] * ( vel_trm_TE_BW + vel );
            pBW[x] = omega_trm * pBW[x] + omega_w2 * rho[x] * ( vel_trm_TE_BW - vel );
         )

         X_LOOP
         (
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTS[x] = omega_trm * pTS[x] + omega_w2 * rho[x] * ( vel_trm_TS_BN - vel );
            pBN[x] = omega_trm * pBN[x] + omega_w2 * rho[x] * ( vel_trm_TS_BN + vel );
         )

         X_LOOP
         (
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            pTN[x] = omega_trm * pTN[x] + omega_w2 * rho[x] * ( vel_trm_TN_BS + vel );
            pBS[x] = omega_trm * pBS[x] + omega_w2 * rho[x] * ( vel_trm_TN_BS - vel );
         )

         X_LOOP
         (
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            pN[x] = omega_trm * pN[x] + omega_w1 * rho[x] * ( vel_trm_N_S + velY[x] );
            pS[x] = omega_trm * pS[x] + omega_w1 * rho[x] * ( vel_trm_N_S - velY[x] );
         )

         X_LOOP
         (
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            pE[x] = omega_trm * pE[x] + omega_w1 * rho[x] * ( vel_trm_E_W + velX[x] );
            pW[x] = omega_trm * pW[x] + omega_w1 * rho[x] * ( vel_trm_E_W - velX[x] );
         )

         X_LOOP
         (
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            pT[x] = omega_trm * pT[x] + omega_w1 * rho[x] * ( vel_trm_T_B + velZ[x] );
            pB[x] = omega_trm * pB[x] + omega_w1 * rho[x] * ( vel_trm_T_B - velZ[x] );
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t dd_tmp_C  = src->get( x, y, z, Stencil::idx[C]  );
            const real_t dd_tmp_N  = src->get( x, y, z, Stencil::idx[N]  );
            const real_t dd_tmp_S  = src->get( x, y, z, Stencil::idx[S]  );
            const real_t dd_tmp_W  = src->get( x, y, z, Stencil::idx[W]  );
            const real_t dd_tmp_E  = src->get( x, y, z, Stencil::idx[E]  );
            const real_t dd_tmp_T  = src->get( x, y, z, Stencil::idx[T]  );
            const real_t dd_tmp_B  = src->get( x, y, z, Stencil::idx[B]  );
            const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW] );
            const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE] );
            const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW] );
            const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE] );
            const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN] );
            const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS] );
            const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW] );
            const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE] );
            const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN] );
            const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS] );
            const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW] );
            const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE] );

            const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
            const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
            const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

            rho[x] = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;
            const real_t rho_inv = real_t(1) / rho[x];

            velX[x] = rho_inv * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
            velY[x] = rho_inv * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
            velZ[x] = rho_inv * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

            dir_indep_trm[x] = one_third - real_c(0.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

            src->get(x,y,z,Stencil::idx[C]) = omega_trm * dd_tmp_C + omega_w0 * rho[x] * dir_indep_trm[x];
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velY[x];
            const real_t vel_trm_NW_SE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[NW]) = omega_trm * src->get(x,y,z,Stencil::idx[NW]) + omega_w2 * rho[x] * ( vel_trm_NW_SE - vel );
            src->get(x,y,z,Stencil::idx[SE]) = omega_trm * src->get(x,y,z,Stencil::idx[SE]) + omega_w2 * rho[x] * ( vel_trm_NW_SE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velY[x];
            const real_t vel_trm_NE_SW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[NE]) = omega_trm * src->get(x,y,z,Stencil::idx[NE]) + omega_w2 * rho[x] * ( vel_trm_NE_SW + vel );
            src->get(x,y,z,Stencil::idx[SW]) = omega_trm * src->get(x,y,z,Stencil::idx[SW]) + omega_w2 * rho[x] * ( vel_trm_NE_SW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] - velZ[x];
            const real_t vel_trm_TW_BE = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TW]) = omega_trm * src->get(x,y,z,Stencil::idx[TW]) + omega_w2 * rho[x] * ( vel_trm_TW_BE - vel );
            src->get(x,y,z,Stencil::idx[BE]) = omega_trm * src->get(x,y,z,Stencil::idx[BE]) + omega_w2 * rho[x] * ( vel_trm_TW_BE + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velX[x] + velZ[x];
            const real_t vel_trm_TE_BW = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TE]) = omega_trm * src->get(x,y,z,Stencil::idx[TE]) + omega_w2 * rho[x] * ( vel_trm_TE_BW + vel );
            src->get(x,y,z,Stencil::idx[BW]) = omega_trm * src->get(x,y,z,Stencil::idx[BW]) + omega_w2 * rho[x] * ( vel_trm_TE_BW - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] - velZ[x];
            const real_t vel_trm_TS_BN = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TS]) = omega_trm * src->get(x,y,z,Stencil::idx[TS]) + omega_w2 * rho[x] * ( vel_trm_TS_BN - vel );
            src->get(x,y,z,Stencil::idx[BN]) = omega_trm * src->get(x,y,z,Stencil::idx[BN]) + omega_w2 * rho[x] * ( vel_trm_TS_BN + vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel = velY[x] + velZ[x];
            const real_t vel_trm_TN_BS = dir_indep_trm[x] + real_c(1.5) * vel * vel;

            src->get(x,y,z,Stencil::idx[TN]) = omega_trm * src->get(x,y,z,Stencil::idx[TN]) + omega_w2 * rho[x] * ( vel_trm_TN_BS + vel );
            src->get(x,y,z,Stencil::idx[BS]) = omega_trm * src->get(x,y,z,Stencil::idx[BS]) + omega_w2 * rho[x] * ( vel_trm_TN_BS - vel );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_N_S = dir_indep_trm[x] + real_c(1.5) * velY[x] * velY[x];

            src->get(x,y,z,Stencil::idx[N]) = omega_trm * src->get(x,y,z,Stencil::idx[N]) + omega_w1 * rho[x] * ( vel_trm_N_S + velY[x] );
            src->get(x,y,z,Stencil::idx[S]) = omega_trm * src->get(x,y,z,Stencil::idx[S]) + omega_w1 * rho[x] * ( vel_trm_N_S - velY[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_E_W = dir_indep_trm[x] + real_c(1.5) * velX[x] * velX[x];

            src->get(x,y,z,Stencil::idx[E]) = omega_trm * src->get(x,y,z,Stencil::idx[E]) + omega_w1 * rho[x] * ( vel_trm_E_W + velX[x] );
            src->get(x,y,z,Stencil::idx[W]) = omega_trm * src->get(x,y,z,Stencil::idx[W]) + omega_w1 * rho[x] * ( vel_trm_E_W - velX[x] );
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            const real_t vel_trm_T_B = dir_indep_trm[x] + real_c(1.5) * velZ[x] * velZ[x];

            src->get(x,y,z,Stencil::idx[T]) = omega_trm * src->get(x,y,z,Stencil::idx[T]) + omega_w1 * rho[x] * ( vel_trm_T_B + velZ[x] );
            src->get(x,y,z,Stencil::idx[B]) = omega_trm * src->get(x,y,z,Stencil::idx[B]) + omega_w1 * rho[x] * ( vel_trm_T_B - velZ[x] );
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] rho;
   delete[] dir_indep_trm;

#ifdef _OPENMP
   }
#endif
}



} // namespace lbm
} // namespace walberla


#else // #ifdef __IBMCPP__

#include "lbm/srt/bluegeneq/SplitPureSweep.impl.h"

#endif
