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
//! \file SplitSweep.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/IntelCompilerOptimization.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/LatticeModelBase.h"
#include "lbm/sweeps/FlagFieldSweepBase.h"
#include "lbm/sweeps/Streaming.h"

#include "field/iterators/IteratorMacros.h"

#include <type_traits>



namespace walberla {
namespace lbm {


///////////////////////////////////////////////////////
// Available TRT implementations:                    //
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

template< typename LatticeModel_T, typename FlagField_T >
class SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                        std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                        ! LatticeModel_T::compressible &&
                                                                        std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                        >::type > :
   public FlagFieldSweepBase< LatticeModel_T, FlagField_T >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value), "Only works with TRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename FlagFieldSweepBase<LatticeModel_T,FlagField_T>::PdfField_T  PdfField_T;
   typedef typename LatticeModel_T::Stencil                                     Stencil;

   // block has NO dst pdf field, lbm mask consists of multiple flags
   SplitSweep( const BlockDataID & pdfField, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( pdfField, flagField, lbmMask ) {}

   // every block has a dedicated dst pdf field, lbm mask consists of multiple flags
   SplitSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( src, dst, flagField, lbmMask ) {}

   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                       ! LatticeModel_T::compressible &&
                                                                       std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                       >::type
   >::operator()( IBlock * const block )                                               
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), 1 );

   // constants used during stream/collide

   const real_t lambda_e =  src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d =  src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT feq_common = new real_t[ uint_c(xSize) ];

   bool * WALBERLA_RESTRICT perform_lbm = new bool[ uint_c(xSize) ];

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
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
               const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
               const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

               const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;

               velX[x] = velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x];
               velY[x] = velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x];
               velZ[x] = velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x];

               feq_common[x] = rho - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               dC[x] = pC[x] * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         )

         real_t * WALBERLA_RESTRICT dNE = &dst->get(0,y,z,Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT dSW = &dst->get(0,y,z,Stencil::idx[SW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( pNE[x] + pSW[x] - fac2 * velXPY * velXPY - t2x2 * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( pNE[x] - pSW[x] - real_t(3.0) * t2x2 * velXPY );

               dNE[x] = pNE[x] - sym_NE_SW - asym_NE_SW;
               dSW[x] = pSW[x] - sym_NE_SW + asym_NE_SW;
            }
         )

         real_t * WALBERLA_RESTRICT dSE = &dst->get(0,y,z,Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT dNW = &dst->get(0,y,z,Stencil::idx[NW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( pSE[x] + pNW[x] - fac2 * velXMY * velXMY - t2x2 * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( pSE[x] - pNW[x] - real_t(3.0) * t2x2 * velXMY );

               dSE[x] = pSE[x] - sym_SE_NW - asym_SE_NW;
               dNW[x] = pNW[x] - sym_SE_NW + asym_SE_NW;
            }
         )

         real_t * WALBERLA_RESTRICT dTE = &dst->get(0,y,z,Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT dBW = &dst->get(0,y,z,Stencil::idx[BW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( pTE[x] + pBW[x] - fac2 * velXPZ * velXPZ - t2x2 * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( pTE[x] - pBW[x] - real_t(3.0) * t2x2 * velXPZ );

               dTE[x] = pTE[x] - sym_TE_BW - asym_TE_BW;
               dBW[x] = pBW[x] - sym_TE_BW + asym_TE_BW;
            }
         )

         real_t * WALBERLA_RESTRICT dBE = &dst->get(0,y,z,Stencil::idx[BE]);
         real_t * WALBERLA_RESTRICT dTW = &dst->get(0,y,z,Stencil::idx[TW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( pBE[x] + pTW[x] - fac2 * velXMZ * velXMZ - t2x2 * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( pBE[x] - pTW[x] - real_t(3.0) * t2x2 * velXMZ );

               dBE[x] = pBE[x] - sym_BE_TW - asym_BE_TW;
               dTW[x] = pTW[x] - sym_BE_TW + asym_BE_TW;
            }
         )

         real_t * WALBERLA_RESTRICT dTN = &dst->get(0,y,z,Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT dBS = &dst->get(0,y,z,Stencil::idx[BS]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( pTN[x] + pBS[x] - fac2 * velYPZ * velYPZ - t2x2 * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( pTN[x] - pBS[x] - real_t(3.0) * t2x2 * velYPZ );

               dTN[x] = pTN[x] - sym_TN_BS - asym_TN_BS;
               dBS[x] = pBS[x] - sym_TN_BS + asym_TN_BS;
            }
         )

         real_t * WALBERLA_RESTRICT dBN = &dst->get(0,y,z,Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT dTS = &dst->get(0,y,z,Stencil::idx[TS]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( pBN[x] + pTS[x] - fac2 * velYMZ * velYMZ - t2x2 * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( pBN[x] - pTS[x] - real_t(3.0) * t2x2 * velYMZ );

               dBN[x] = pBN[x] - sym_BN_TS - asym_BN_TS;
               dTS[x] = pTS[x] - sym_BN_TS + asym_BN_TS;
            }
         )

         real_t * WALBERLA_RESTRICT dN = &dst->get(0,y,z,Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT dS = &dst->get(0,y,z,Stencil::idx[S]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_N_S = lambda_e_scaled * ( pN[x] + pS[x] - fac1 * velY[x] * velY[x] - t1x2 * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( pN[x] - pS[x] - real_t(3.0) * t1x2 * velY[x] );

               dN[x] = pN[x] - sym_N_S - asym_N_S;
               dS[x] = pS[x] - sym_N_S + asym_N_S;
            }
         )

         real_t * WALBERLA_RESTRICT dE = &dst->get(0,y,z,Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT dW = &dst->get(0,y,z,Stencil::idx[W]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_E_W = lambda_e_scaled * ( pE[x] + pW[x] - fac1 * velX[x] * velX[x] - t1x2 * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( pE[x] - pW[x] - real_t(3.0) * t1x2 * velX[x] );

               dE[x] = pE[x] - sym_E_W - asym_E_W;
               dW[x] = pW[x] - sym_E_W + asym_E_W;
            }
         )

         real_t * WALBERLA_RESTRICT dT = &dst->get(0,y,z,Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT dB = &dst->get(0,y,z,Stencil::idx[B]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_T_B = lambda_e_scaled * ( pT[x] + pB[x] - fac1 * velZ[x] * velZ[x] - t1x2 * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( pT[x] - pB[x] - real_t(3.0) * t1x2 * velZ[x] );

               dT[x] = pT[x] - sym_T_B - asym_T_B;
               dB[x] = pB[x] - sym_T_B + asym_T_B;
            }
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
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

               feq_common[x] = rho - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               dst->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_NE = src->get(x-1, y-1, z, Stencil::idx[NE]);
               const real_t dd_tmp_SW = src->get(x+1, y+1, z, Stencil::idx[SW]);

               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( dd_tmp_NE + dd_tmp_SW - fac2 * velXPY * velXPY - t2x2 * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2 * velXPY );

               dst->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
               dst->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_SE = src->get(x-1, y+1, z, Stencil::idx[SE]);
               const real_t dd_tmp_NW = src->get(x+1, y-1, z, Stencil::idx[NW]);

               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( dd_tmp_SE + dd_tmp_NW - fac2 * velXMY * velXMY - t2x2 * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2 * velXMY );

               dst->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
               dst->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TE = src->get(x-1, y, z-1, Stencil::idx[TE]);
               const real_t dd_tmp_BW = src->get(x+1, y, z+1, Stencil::idx[BW]);

               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( dd_tmp_TE + dd_tmp_BW - fac2 * velXPZ * velXPZ - t2x2 * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2 * velXPZ );

               dst->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
               dst->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BE = src->get(x-1, y, z+1, Stencil::idx[BE]);
               const real_t dd_tmp_TW = src->get(x+1, y, z-1, Stencil::idx[TW]);

               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( dd_tmp_BE + dd_tmp_TW - fac2 * velXMZ * velXMZ - t2x2 * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2 * velXMZ );

               dst->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
               dst->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TN = src->get(x, y-1, z-1, Stencil::idx[TN]);
               const real_t dd_tmp_BS = src->get(x, y+1, z+1, Stencil::idx[BS]);

               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( dd_tmp_TN + dd_tmp_BS - fac2 * velYPZ * velYPZ - t2x2 * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2 * velYPZ );

               dst->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
               dst->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BN = src->get(x, y-1, z+1, Stencil::idx[BN]);
               const real_t dd_tmp_TS = src->get(x, y+1, z-1, Stencil::idx[TS]);

               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( dd_tmp_BN + dd_tmp_TS - fac2 * velYMZ * velYMZ - t2x2 * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2 * velYMZ );

               dst->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
               dst->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_N  = src->get(x, y-1, z, Stencil::idx[N]);
               const real_t dd_tmp_S  = src->get(x, y+1, z, Stencil::idx[S]);

               const real_t  sym_N_S = lambda_e_scaled * ( dd_tmp_N + dd_tmp_S - fac1 * velY[x] * velY[x] - t1x2 * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2 * velY[x] );

               dst->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
               dst->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_E  = src->get(x-1, y, z, Stencil::idx[E]);
               const real_t dd_tmp_W  = src->get(x+1, y, z, Stencil::idx[W]);

               const real_t  sym_E_W = lambda_e_scaled * ( dd_tmp_E + dd_tmp_W - fac1 * velX[x] * velX[x] - t1x2 * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2 * velX[x] );

               dst->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
               dst->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_T  = src->get(x, y, z-1, Stencil::idx[T]);
               const real_t dd_tmp_B  = src->get(x, y, z+1, Stencil::idx[B]);

               const real_t  sym_T_B = lambda_e_scaled * ( dd_tmp_T + dd_tmp_B - fac1 * velZ[x] * velZ[x] - t1x2 * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2 * velZ[x] );

               dst->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
               dst->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
            }
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] feq_common;
   delete[] perform_lbm;

#ifdef _OPENMP
   }
#endif

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                       ! LatticeModel_T::compressible &&
                                                                       std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                       >::type
   >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   Stream< LatticeModel_T, FlagField_T >::execute( src, dst, flagField, lbm, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
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

   PdfField_T * src( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   // constants used during stream/collide

   const real_t lambda_e =  src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d =  src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT feq_common = new real_t[ uint_c(xSize) ];

   bool * WALBERLA_RESTRICT perform_lbm = new bool[ uint_c(xSize) ];

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
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
               const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
               const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

               const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;

               velX[x] = velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x];
               velY[x] = velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x];
               velZ[x] = velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x];

               feq_common[x] = rho - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               pC[x] = pC[x] * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( pNE[x] + pSW[x] - fac2 * velXPY * velXPY - t2x2 * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( pNE[x] - pSW[x] - real_t(3.0) * t2x2 * velXPY );

               pNE[x] = pNE[x] - sym_NE_SW - asym_NE_SW;
               pSW[x] = pSW[x] - sym_NE_SW + asym_NE_SW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( pSE[x] + pNW[x] - fac2 * velXMY * velXMY - t2x2 * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( pSE[x] - pNW[x] - real_t(3.0) * t2x2 * velXMY );

               pSE[x] = pSE[x] - sym_SE_NW - asym_SE_NW;
               pNW[x] = pNW[x] - sym_SE_NW + asym_SE_NW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( pTE[x] + pBW[x] - fac2 * velXPZ * velXPZ - t2x2 * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( pTE[x] - pBW[x] - real_t(3.0) * t2x2 * velXPZ );

               pTE[x] = pTE[x] - sym_TE_BW - asym_TE_BW;
               pBW[x] = pBW[x] - sym_TE_BW + asym_TE_BW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( pBE[x] + pTW[x] - fac2 * velXMZ * velXMZ - t2x2 * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( pBE[x] - pTW[x] - real_t(3.0) * t2x2 * velXMZ );

               pBE[x] = pBE[x] - sym_BE_TW - asym_BE_TW;
               pTW[x] = pTW[x] - sym_BE_TW + asym_BE_TW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( pTN[x] + pBS[x] - fac2 * velYPZ * velYPZ - t2x2 * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( pTN[x] - pBS[x] - real_t(3.0) * t2x2 * velYPZ );

               pTN[x] = pTN[x] - sym_TN_BS - asym_TN_BS;
               pBS[x] = pBS[x] - sym_TN_BS + asym_TN_BS;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( pBN[x] + pTS[x] - fac2 * velYMZ * velYMZ - t2x2 * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( pBN[x] - pTS[x] - real_t(3.0) * t2x2 * velYMZ );

               pBN[x] = pBN[x] - sym_BN_TS - asym_BN_TS;
               pTS[x] = pTS[x] - sym_BN_TS + asym_BN_TS;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_N_S = lambda_e_scaled * ( pN[x] + pS[x] - fac1 * velY[x] * velY[x] - t1x2 * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( pN[x] - pS[x] - real_t(3.0) * t1x2 * velY[x] );

               pN[x] = pN[x] - sym_N_S - asym_N_S;
               pS[x] = pS[x] - sym_N_S + asym_N_S;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_E_W = lambda_e_scaled * ( pE[x] + pW[x] - fac1 * velX[x] * velX[x] - t1x2 * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( pE[x] - pW[x] - real_t(3.0) * t1x2 * velX[x] );

               pE[x] = pE[x] - sym_E_W - asym_E_W;
               pW[x] = pW[x] - sym_E_W + asym_E_W;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_T_B = lambda_e_scaled * ( pT[x] + pB[x] - fac1 * velZ[x] * velZ[x] - t1x2 * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( pT[x] - pB[x] - real_t(3.0) * t1x2 * velZ[x] );

               pT[x] = pT[x] - sym_T_B - asym_T_B;
               pB[x] = pB[x] - sym_T_B + asym_T_B;
            }
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
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

               feq_common[x] = rho - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               src->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE]);
               const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW]);

               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( dd_tmp_NE + dd_tmp_SW - fac2 * velXPY * velXPY - t2x2 * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2 * velXPY );

               src->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
               src->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE]);
               const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW]);

               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( dd_tmp_SE + dd_tmp_NW - fac2 * velXMY * velXMY - t2x2 * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2 * velXMY );

               src->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
               src->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE]);
               const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW]);

               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( dd_tmp_TE + dd_tmp_BW - fac2 * velXPZ * velXPZ - t2x2 * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2 * velXPZ );

               src->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
               src->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE]);
               const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW]);

               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( dd_tmp_BE + dd_tmp_TW - fac2 * velXMZ * velXMZ - t2x2 * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2 * velXMZ );

               src->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
               src->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN]);
               const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS]);

               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( dd_tmp_TN + dd_tmp_BS - fac2 * velYPZ * velYPZ - t2x2 * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2 * velYPZ );

               src->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
               src->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN]);
               const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS]);

               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( dd_tmp_BN + dd_tmp_TS - fac2 * velYMZ * velYMZ - t2x2 * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2 * velYMZ );

               src->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
               src->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_N  = src->get( x, y, z, Stencil::idx[N]);
               const real_t dd_tmp_S  = src->get( x, y, z, Stencil::idx[S]);

               const real_t  sym_N_S = lambda_e_scaled * ( dd_tmp_N + dd_tmp_S - fac1 * velY[x] * velY[x] - t1x2 * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2 * velY[x] );

               src->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
               src->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_E  = src->get( x, y, z, Stencil::idx[E]);
               const real_t dd_tmp_W  = src->get( x, y, z, Stencil::idx[W]);

               const real_t  sym_E_W = lambda_e_scaled * ( dd_tmp_E + dd_tmp_W - fac1 * velX[x] * velX[x] - t1x2 * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2 * velX[x] );

               src->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
               src->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_T  = src->get( x, y, z, Stencil::idx[T]);
               const real_t dd_tmp_B  = src->get( x, y, z, Stencil::idx[B]);

               const real_t  sym_T_B = lambda_e_scaled * ( dd_tmp_T + dd_tmp_B - fac1 * velZ[x] * velZ[x] - t1x2 * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2 * velZ[x] );

               src->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
               src->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
            }
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] feq_common;
   delete[] perform_lbm;

#ifdef _OPENMP
   }
#endif
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T, typename FlagField_T >
class SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                        std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                        LatticeModel_T::compressible &&
                                                                        std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                        >::type > :
   public FlagFieldSweepBase< LatticeModel_T, FlagField_T >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value), "Only works with TRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename FlagFieldSweepBase<LatticeModel_T,FlagField_T>::PdfField_T  PdfField_T;
   typedef typename LatticeModel_T::Stencil                                     Stencil;

   // block has NO dst pdf field, lbm mask consists of multiple flags
   SplitSweep( const BlockDataID & pdfField, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( pdfField, flagField, lbmMask ) {}

   // every block has a dedicated dst pdf field, lbm mask consists of multiple flags
   SplitSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( src, dst, flagField, lbmMask ) {}

   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                       LatticeModel_T::compressible &&
                                                                       std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                       >::type
   >::operator()( IBlock * const block )
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), 1 );

   // constants used during stream/collide

   const real_t lambda_e =  src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d =  src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2_0 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT t1x2 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT t2x2 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT fac1 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT fac2 = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT feq_common = new real_t[ uint_c(xSize) ];

   bool * WALBERLA_RESTRICT perform_lbm = new bool[ uint_c(xSize) ];

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
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
               const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
               const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

               const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;
               const real_t invRho = real_t(1.0) / rho;

               velX[x] = invRho * ( velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x] );
               velY[x] = invRho * ( velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x] );
               velZ[x] = invRho * ( velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x] );

               t1x2[x] = t1x2_0 * rho;
               t2x2[x] = t2x2_0 * rho;
               fac1[x] = t1x2_0 * rho * inv2csq2;
               fac2[x] = t2x2_0 * rho * inv2csq2;

               feq_common[x] = real_t(1.0) - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               dC[x] = pC[x] * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         )

         real_t * WALBERLA_RESTRICT dNE = &dst->get(0,y,z,Stencil::idx[NE]);
         real_t * WALBERLA_RESTRICT dSW = &dst->get(0,y,z,Stencil::idx[SW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( pNE[x] + pSW[x] - fac2[x] * velXPY * velXPY - t2x2[x] * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( pNE[x] - pSW[x] - real_t(3.0) * t2x2[x] * velXPY );

               dNE[x] = pNE[x] - sym_NE_SW - asym_NE_SW;
               dSW[x] = pSW[x] - sym_NE_SW + asym_NE_SW;
            }
         )

         real_t * WALBERLA_RESTRICT dSE = &dst->get(0,y,z,Stencil::idx[SE]);
         real_t * WALBERLA_RESTRICT dNW = &dst->get(0,y,z,Stencil::idx[NW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( pSE[x] + pNW[x] - fac2[x] * velXMY * velXMY - t2x2[x] * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( pSE[x] - pNW[x] - real_t(3.0) * t2x2[x] * velXMY );

               dSE[x] = pSE[x] - sym_SE_NW - asym_SE_NW;
               dNW[x] = pNW[x] - sym_SE_NW + asym_SE_NW;
            }
         )

         real_t * WALBERLA_RESTRICT dTE = &dst->get(0,y,z,Stencil::idx[TE]);
         real_t * WALBERLA_RESTRICT dBW = &dst->get(0,y,z,Stencil::idx[BW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( pTE[x] + pBW[x] - fac2[x] * velXPZ * velXPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( pTE[x] - pBW[x] - real_t(3.0) * t2x2[x] * velXPZ );

               dTE[x] = pTE[x] - sym_TE_BW - asym_TE_BW;
               dBW[x] = pBW[x] - sym_TE_BW + asym_TE_BW;
            }
         )

         real_t * WALBERLA_RESTRICT dBE = &dst->get(0,y,z,Stencil::idx[BE]);
         real_t * WALBERLA_RESTRICT dTW = &dst->get(0,y,z,Stencil::idx[TW]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( pBE[x] + pTW[x] - fac2[x] * velXMZ * velXMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( pBE[x] - pTW[x] - real_t(3.0) * t2x2[x] * velXMZ );

               dBE[x] = pBE[x] - sym_BE_TW - asym_BE_TW;
               dTW[x] = pTW[x] - sym_BE_TW + asym_BE_TW;
            }
         )

         real_t * WALBERLA_RESTRICT dTN = &dst->get(0,y,z,Stencil::idx[TN]);
         real_t * WALBERLA_RESTRICT dBS = &dst->get(0,y,z,Stencil::idx[BS]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( pTN[x] + pBS[x] - fac2[x] * velYPZ * velYPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( pTN[x] - pBS[x] - real_t(3.0) * t2x2[x] * velYPZ );

               dTN[x] = pTN[x] - sym_TN_BS - asym_TN_BS;
               dBS[x] = pBS[x] - sym_TN_BS + asym_TN_BS;
            }
         )

         real_t * WALBERLA_RESTRICT dBN = &dst->get(0,y,z,Stencil::idx[BN]);
         real_t * WALBERLA_RESTRICT dTS = &dst->get(0,y,z,Stencil::idx[TS]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( pBN[x] + pTS[x] - fac2[x] * velYMZ * velYMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( pBN[x] - pTS[x] - real_t(3.0) * t2x2[x] * velYMZ );

               dBN[x] = pBN[x] - sym_BN_TS - asym_BN_TS;
               dTS[x] = pTS[x] - sym_BN_TS + asym_BN_TS;
            }
         )

         real_t * WALBERLA_RESTRICT dN = &dst->get(0,y,z,Stencil::idx[N]);
         real_t * WALBERLA_RESTRICT dS = &dst->get(0,y,z,Stencil::idx[S]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_N_S = lambda_e_scaled * ( pN[x] + pS[x] - fac1[x] * velY[x] * velY[x] - t1x2[x] * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( pN[x] - pS[x] - real_t(3.0) * t1x2[x] * velY[x] );

               dN[x] = pN[x] - sym_N_S - asym_N_S;
               dS[x] = pS[x] - sym_N_S + asym_N_S;
            }
         )

         real_t * WALBERLA_RESTRICT dE = &dst->get(0,y,z,Stencil::idx[E]);
         real_t * WALBERLA_RESTRICT dW = &dst->get(0,y,z,Stencil::idx[W]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_E_W = lambda_e_scaled * ( pE[x] + pW[x] - fac1[x] * velX[x] * velX[x] - t1x2[x] * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( pE[x] - pW[x] - real_t(3.0) * t1x2[x] * velX[x] );

               dE[x] = pE[x] - sym_E_W - asym_E_W;
               dW[x] = pW[x] - sym_E_W + asym_E_W;
            }
         )

         real_t * WALBERLA_RESTRICT dT = &dst->get(0,y,z,Stencil::idx[T]);
         real_t * WALBERLA_RESTRICT dB = &dst->get(0,y,z,Stencil::idx[B]);

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_T_B = lambda_e_scaled * ( pT[x] + pB[x] - fac1[x] * velZ[x] * velZ[x] - t1x2[x] * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( pT[x] - pB[x] - real_t(3.0) * t1x2[x] * velZ[x] );

               dT[x] = pT[x] - sym_T_B - asym_T_B;
               dB[x] = pB[x] - sym_T_B + asym_T_B;
            }
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
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
               const real_t invRho = real_t(1.0) / rho;

               velX[x] = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
               velY[x] = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
               velZ[x] = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

               t1x2[x] = t1x2_0 * rho;
               t2x2[x] = t2x2_0 * rho;
               fac1[x] = t1x2_0 * rho * inv2csq2;
               fac2[x] = t2x2_0 * rho * inv2csq2;

               feq_common[x] = real_t(1.0) - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               dst->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_NE = src->get(x-1, y-1, z, Stencil::idx[NE]);
               const real_t dd_tmp_SW = src->get(x+1, y+1, z, Stencil::idx[SW]);

               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( dd_tmp_NE + dd_tmp_SW - fac2[x] * velXPY * velXPY - t2x2[x] * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2[x] * velXPY );

               dst->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
               dst->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_SE = src->get(x-1, y+1, z, Stencil::idx[SE]);
               const real_t dd_tmp_NW = src->get(x+1, y-1, z, Stencil::idx[NW]);

               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( dd_tmp_SE + dd_tmp_NW - fac2[x] * velXMY * velXMY - t2x2[x] * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2[x] * velXMY );

               dst->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
               dst->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TE = src->get(x-1, y, z-1, Stencil::idx[TE]);
               const real_t dd_tmp_BW = src->get(x+1, y, z+1, Stencil::idx[BW]);

               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( dd_tmp_TE + dd_tmp_BW - fac2[x] * velXPZ * velXPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2[x] * velXPZ );

               dst->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
               dst->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BE = src->get(x-1, y, z+1, Stencil::idx[BE]);
               const real_t dd_tmp_TW = src->get(x+1, y, z-1, Stencil::idx[TW]);

               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( dd_tmp_BE + dd_tmp_TW - fac2[x] * velXMZ * velXMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2[x] * velXMZ );

               dst->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
               dst->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TN = src->get(x, y-1, z-1, Stencil::idx[TN]);
               const real_t dd_tmp_BS = src->get(x, y+1, z+1, Stencil::idx[BS]);

               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( dd_tmp_TN + dd_tmp_BS - fac2[x] * velYPZ * velYPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2[x] * velYPZ );

               dst->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
               dst->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BN = src->get(x, y-1, z+1, Stencil::idx[BN]);
               const real_t dd_tmp_TS = src->get(x, y+1, z-1, Stencil::idx[TS]);

               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( dd_tmp_BN + dd_tmp_TS - fac2[x] * velYMZ * velYMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2[x] * velYMZ );

               dst->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
               dst->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_N  = src->get(x, y-1, z, Stencil::idx[N]);
               const real_t dd_tmp_S  = src->get(x, y+1, z, Stencil::idx[S]);

               const real_t  sym_N_S = lambda_e_scaled * ( dd_tmp_N + dd_tmp_S - fac1[x] * velY[x] * velY[x] - t1x2[x] * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2[x] * velY[x] );

               dst->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
               dst->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_E  = src->get(x-1, y, z, Stencil::idx[E]);
               const real_t dd_tmp_W  = src->get(x+1, y, z, Stencil::idx[W]);

               const real_t  sym_E_W = lambda_e_scaled * ( dd_tmp_E + dd_tmp_W - fac1[x] * velX[x] * velX[x] - t1x2[x] * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2[x] * velX[x] );

               dst->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
               dst->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_T  = src->get(x, y, z-1, Stencil::idx[T]);
               const real_t dd_tmp_B  = src->get(x, y, z+1, Stencil::idx[B]);

               const real_t  sym_T_B = lambda_e_scaled * ( dd_tmp_T + dd_tmp_B - fac1[x] * velZ[x] * velZ[x] - t1x2[x] * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2[x] * velZ[x] );

               dst->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
               dst->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
            }
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] t1x2;
   delete[] t2x2;
   delete[] fac1;
   delete[] fac2;
   delete[] feq_common;
   delete[] perform_lbm;

#ifdef _OPENMP
   }
#endif

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                       LatticeModel_T::compressible &&
                                                                       std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                       >::type
   >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   Stream< LatticeModel_T, FlagField_T >::execute( src, dst, flagField, lbm, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T, typename FlagField_T >
void SplitSweep< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
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

   PdfField_T * src( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   // constants used during stream/collide

   const real_t lambda_e =  src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d =  src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2_0 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   // loop constants

   const cell_idx_t xSize = cell_idx_c( src->xSize() );

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif
   // temporaries, calculated by the first innermost loop

   real_t * WALBERLA_RESTRICT velX = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velY = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT velZ = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT t1x2 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT t2x2 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT fac1 = new real_t[ uint_c(xSize) ];
   real_t * WALBERLA_RESTRICT fac2 = new real_t[ uint_c(xSize) ];

   real_t * WALBERLA_RESTRICT feq_common = new real_t[ uint_c(xSize) ];

   bool * WALBERLA_RESTRICT perform_lbm = new bool[ uint_c(xSize) ];

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
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const real_t velX_trm = pE[x] + pNE[x] + pSE[x] + pTE[x] + pBE[x];
               const real_t velY_trm = pN[x] + pNW[x] + pTN[x] + pBN[x];
               const real_t velZ_trm = pT[x] + pTS[x] + pTW[x];

               const real_t rho = pC[x] + pS[x] + pW[x] + pB[x] + pSW[x] + pBS[x] + pBW[x] + velX_trm + velY_trm + velZ_trm;
               const real_t invRho = real_t(1.0) / rho;

               velX[x] = invRho * ( velX_trm - pW[x]  - pNW[x] - pSW[x] - pTW[x] - pBW[x] );
               velY[x] = invRho * ( velY_trm + pNE[x] - pS[x]  - pSW[x] - pSE[x] - pTS[x] - pBS[x] );
               velZ[x] = invRho * ( velZ_trm + pTN[x] + pTE[x] - pB[x]  - pBN[x] - pBS[x] - pBW[x] - pBE[x] );

               t1x2[x] = t1x2_0 * rho;
               t2x2[x] = t2x2_0 * rho;
               fac1[x] = t1x2_0 * rho * inv2csq2;
               fac2[x] = t2x2_0 * rho * inv2csq2;

               feq_common[x] = real_t(1.0) - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               pC[x] = pC[x] * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( pNE[x] + pSW[x] - fac2[x] * velXPY * velXPY - t2x2[x] * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( pNE[x] - pSW[x] - real_t(3.0) * t2x2[x] * velXPY );

               pNE[x] = pNE[x] - sym_NE_SW - asym_NE_SW;
               pSW[x] = pSW[x] - sym_NE_SW + asym_NE_SW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( pSE[x] + pNW[x] - fac2[x] * velXMY * velXMY - t2x2[x] * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( pSE[x] - pNW[x] - real_t(3.0) * t2x2[x] * velXMY );

               pSE[x] = pSE[x] - sym_SE_NW - asym_SE_NW;
               pNW[x] = pNW[x] - sym_SE_NW + asym_SE_NW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( pTE[x] + pBW[x] - fac2[x] * velXPZ * velXPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( pTE[x] - pBW[x] - real_t(3.0) * t2x2[x] * velXPZ );

               pTE[x] = pTE[x] - sym_TE_BW - asym_TE_BW;
               pBW[x] = pBW[x] - sym_TE_BW + asym_TE_BW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( pBE[x] + pTW[x] - fac2[x] * velXMZ * velXMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( pBE[x] - pTW[x] - real_t(3.0) * t2x2[x] * velXMZ );

               pBE[x] = pBE[x] - sym_BE_TW - asym_BE_TW;
               pTW[x] = pTW[x] - sym_BE_TW + asym_BE_TW;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( pTN[x] + pBS[x] - fac2[x] * velYPZ * velYPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( pTN[x] - pBS[x] - real_t(3.0) * t2x2[x] * velYPZ );

               pTN[x] = pTN[x] - sym_TN_BS - asym_TN_BS;
               pBS[x] = pBS[x] - sym_TN_BS + asym_TN_BS;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( pBN[x] + pTS[x] - fac2[x] * velYMZ * velYMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( pBN[x] - pTS[x] - real_t(3.0) * t2x2[x] * velYMZ );

               pBN[x] = pBN[x] - sym_BN_TS - asym_BN_TS;
               pTS[x] = pTS[x] - sym_BN_TS + asym_BN_TS;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_N_S = lambda_e_scaled * ( pN[x] + pS[x] - fac1[x] * velY[x] * velY[x] - t1x2[x] * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( pN[x] - pS[x] - real_t(3.0) * t1x2[x] * velY[x] );

               pN[x] = pN[x] - sym_N_S - asym_N_S;
               pS[x] = pS[x] - sym_N_S + asym_N_S;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_E_W = lambda_e_scaled * ( pE[x] + pW[x] - fac1[x] * velX[x] * velX[x] - t1x2[x] * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( pE[x] - pW[x] - real_t(3.0) * t1x2[x] * velX[x] );

               pE[x] = pE[x] - sym_E_W - asym_E_W;
               pW[x] = pW[x] - sym_E_W + asym_E_W;
            }
         )

         X_LOOP
         (
            if( perform_lbm[x] )
            {
               const real_t  sym_T_B = lambda_e_scaled * ( pT[x] + pB[x] - fac1[x] * velZ[x] * velZ[x] - t1x2[x] * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( pT[x] - pB[x] - real_t(3.0) * t1x2[x] * velZ[x] );

               pT[x] = pT[x] - sym_T_B - asym_T_B;
               pB[x] = pB[x] - sym_T_B + asym_T_B;
            }
         )

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }
   else // ==> src->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_YZ_OMP( src, omp for schedule(static),

         using namespace stencil;

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
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
               const real_t invRho = real_t(1.0) / rho;

               velX[x] = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
               velY[x] = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
               velZ[x] = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

               t1x2[x] = t1x2_0 * rho;
               t2x2[x] = t2x2_0 * rho;
               fac1[x] = t1x2_0 * rho * inv2csq2;
               fac2[x] = t2x2_0 * rho * inv2csq2;

               feq_common[x] = real_t(1.0) - real_t(1.5) * ( velX[x] * velX[x] + velY[x] * velY[x] + velZ[x] * velZ[x] );

               src->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common[x];

               perform_lbm[x] = true;
            }
            else perform_lbm[x] = false;
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE]);
               const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW]);

               const real_t velXPY = velX[x] + velY[x];
               const real_t  sym_NE_SW = lambda_e_scaled * ( dd_tmp_NE + dd_tmp_SW - fac2[x] * velXPY * velXPY - t2x2[x] * feq_common[x] );
               const real_t asym_NE_SW = lambda_d_scaled * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2[x] * velXPY );

               src->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
               src->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE]);
               const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW]);

               const real_t velXMY = velX[x] - velY[x];
               const real_t  sym_SE_NW = lambda_e_scaled * ( dd_tmp_SE + dd_tmp_NW - fac2[x] * velXMY * velXMY - t2x2[x] * feq_common[x] );
               const real_t asym_SE_NW = lambda_d_scaled * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2[x] * velXMY );

               src->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
               src->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE]);
               const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW]);

               const real_t velXPZ = velX[x] + velZ[x];
               const real_t  sym_TE_BW = lambda_e_scaled * ( dd_tmp_TE + dd_tmp_BW - fac2[x] * velXPZ * velXPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TE_BW = lambda_d_scaled * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2[x] * velXPZ );

               src->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
               src->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE]);
               const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW]);

               const real_t velXMZ = velX[x] - velZ[x];
               const real_t  sym_BE_TW = lambda_e_scaled * ( dd_tmp_BE + dd_tmp_TW - fac2[x] * velXMZ * velXMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BE_TW = lambda_d_scaled * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2[x] * velXMZ );

               src->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
               src->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN]);
               const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS]);

               const real_t velYPZ = velY[x] + velZ[x];
               const real_t  sym_TN_BS = lambda_e_scaled * ( dd_tmp_TN + dd_tmp_BS - fac2[x] * velYPZ * velYPZ - t2x2[x] * feq_common[x] );
               const real_t asym_TN_BS = lambda_d_scaled * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2[x] * velYPZ );

               src->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
               src->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN]);
               const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS]);

               const real_t velYMZ = velY[x] - velZ[x];
               const real_t  sym_BN_TS = lambda_e_scaled * ( dd_tmp_BN + dd_tmp_TS - fac2[x] * velYMZ * velYMZ - t2x2[x] * feq_common[x] );
               const real_t asym_BN_TS = lambda_d_scaled * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2[x] * velYMZ );

               src->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
               src->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_N  = src->get( x, y, z, Stencil::idx[N]);
               const real_t dd_tmp_S  = src->get( x, y, z, Stencil::idx[S]);

               const real_t  sym_N_S = lambda_e_scaled * ( dd_tmp_N + dd_tmp_S - fac1[x] * velY[x] * velY[x] - t1x2[x] * feq_common[x] );
               const real_t asym_N_S = lambda_d_scaled * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2[x] * velY[x] );

               src->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
               src->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_E  = src->get( x, y, z, Stencil::idx[E]);
               const real_t dd_tmp_W  = src->get( x, y, z, Stencil::idx[W]);

               const real_t  sym_E_W = lambda_e_scaled * ( dd_tmp_E + dd_tmp_W - fac1[x] * velX[x] * velX[x] - t1x2[x] * feq_common[x] );
               const real_t asym_E_W = lambda_d_scaled * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2[x] * velX[x] );

               src->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
               src->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;
            }
         }

         for( cell_idx_t x = 0; x != xSize; ++x )
         {
            if( perform_lbm[x] )
            {
               const real_t dd_tmp_T  = src->get( x, y, z, Stencil::idx[T]);
               const real_t dd_tmp_B  = src->get( x, y, z, Stencil::idx[B]);

               const real_t  sym_T_B = lambda_e_scaled * ( dd_tmp_T + dd_tmp_B - fac1[x] * velZ[x] * velZ[x] - t1x2[x] * feq_common[x] );
               const real_t asym_T_B = lambda_d_scaled * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2[x] * velZ[x] );

               src->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
               src->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
            }
         }

      ) // WALBERLA_FOR_ALL_CELLS_YZ_OMP
   }

   delete[] velX;
   delete[] velY;
   delete[] velZ;
   delete[] t1x2;
   delete[] t2x2;
   delete[] fac1;
   delete[] fac2;
   delete[] feq_common;
   delete[] perform_lbm;

#ifdef _OPENMP
   }
#endif
}



} // namespace lbm
} // namespace walberla
