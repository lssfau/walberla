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
//! \file CellwiseSweep.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/field/DensityVelocityCallback.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"
#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include <type_traits>


namespace walberla {
namespace lbm {

// For a generic implementation of the TRT sweep jump to the end of this file!



//////////////////////////
// D2Q9 SPECIALIZATIONS //
//////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(4.0) / real_t( 9.0);                // 1/3      for C
   const real_t t1x2 = real_t(1.0) / real_t( 9.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, real_t(0) ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(4.0) / real_t( 9.0);                // 1/3      for C
   const real_t t1x2 = real_t(1.0) / real_t( 9.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, real_t(0) ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7' at the end of this file!






///////////////////////////
// D3Q19 SPECIALIZATIONS //
///////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2_0 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t feq_common = real_t(1.0) - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common;

         const real_t t2x2 = t2x2_0 * rho;
         const real_t fac2 = t2x2 * inv2csq2;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t t1x2 = t1x2_0 * rho;
         const real_t fac1 = t1x2 * inv2csq2;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(1.0) / real_t(3.0);                 // 1/3      for C
   const real_t t1x2_0 = real_t(1.0) / real_t(18.0) * real_t(2.0);  // 1/18 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(36.0) * real_t(2.0);  // 1/36 * 2 else

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t feq_common = real_t(1.0) - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common;

         const real_t t2x2 = t2x2_0 * rho;
         const real_t fac2 = t2x2 * inv2csq2;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t t1x2 = t1x2_0 * rho;
         const real_t fac1 = t1x2 * inv2csq2;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common; // no force term

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);
         
         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW + three_w2 * (  force[0] + force[1] );
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW + three_w2 * ( -force[0] - force[1] );

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW + three_w2 * (  force[0] - force[1] );
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW + three_w2 * (  force[1] - force[0] );

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW + three_w2 * (  force[0] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW + three_w2 * ( -force[0] - force[2] );

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW + three_w2 * (  force[0] - force[2] );
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW + three_w2 * (  force[2] - force[0] );

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS + three_w2 * (  force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS + three_w2 * ( -force[1] - force[2] );

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS + three_w2 * (  force[1] - force[2] );
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS + three_w2 * (  force[2] - force[1] );

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S + three_w1 * force[1];
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S - three_w1 * force[1];

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W + three_w1 * force[0];
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W - three_w1 * force[0];

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B + three_w1 * force[2];
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B - three_w1 * force[2];
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common; // no force term

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);
         
         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW + three_w2 * (  force[0] + force[1] );
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW + three_w2 * ( -force[0] - force[1] );

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW + three_w2 * (  force[0] - force[1] );
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW + three_w2 * (  force[1] - force[0] );

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW + three_w2 * (  force[0] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW + three_w2 * ( -force[0] - force[2] );

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW + three_w2 * (  force[0] - force[2] );
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW + three_w2 * (  force[2] - force[0] );

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS + three_w2 * (  force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS + three_w2 * ( -force[1] - force[2] );

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS + three_w2 * (  force[1] - force[2] );
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS + three_w2 * (  force[2] - force[1] );

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S + three_w1 * force[1];
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S - three_w1 * force[1];

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W + three_w1 * force[0];
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W - three_w1 * force[0];

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B + three_w1 * force[2];
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B - three_w1 * force[2];
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3' at the end of this file!






///////////////////////////
// D3Q27 SPECIALIZATIONS //
///////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation>::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2 = real_t(2.0) / real_t(27.0) * real_t(2.0);  // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(54.0) * real_t(2.0);  // 1/54 * 2 else
   const real_t t3x2 = real_t(1.0) / real_t(216.0) * real_t(2.0); // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;
   const real_t fac3     = t3x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations


   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         dst->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW;
         dst->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW;

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         dst->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE;
         dst->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE;

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         dst->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW;
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW;

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         dst->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE;
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE;

      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2 = real_t(2.0) / real_t(27.0) * real_t(2.0);  // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(54.0) * real_t(2.0);  // 1/54 * 2 else
   const real_t t3x2 = real_t(1.0) / real_t(216.0)* real_t(2.0);  // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;
   const real_t fac3     = t3x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         src->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW;
         src->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW;

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         src->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE;
         src->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE;

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         src->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW;
         src->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW;

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         src->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE;
         src->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2_0 = real_t(2.0) / real_t(27.0)  * real_t(2.0); // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(54.0)  * real_t(2.0); // 1/54 * 2 else
   const real_t t3x2_0 = real_t(1.0) / real_t(216.0) * real_t(2.0); // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t feq_common = real_t(1.0) - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common;

         const real_t t2x2 = t2x2_0 * rho;
         const real_t fac2 = t2x2 * inv2csq2;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t t1x2 = t1x2_0 * rho;
         const real_t fac1 = t1x2 * inv2csq2;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;

         const real_t t3x2 = t3x2_0 * rho;
         const real_t fac3 = t3x2 * inv2csq2;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         dst->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW;
         dst->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW;

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         dst->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE;
         dst->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE;

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         dst->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW;
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW;

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         dst->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE;
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5 )
{
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   // common prefactors for calculating the equilibrium parts
   const real_t t0_0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2_0 = real_t(2.0) / real_t(27.0) * real_t(2.0);  // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2_0 = real_t(1.0) / real_t(54.0) * real_t(2.0);  // 1/54 * 2 else
   const real_t t3x2_0 = real_t(1.0) / real_t(216.0) * real_t(2.0); // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t feq_common = real_t(1.0) - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0_0 * rho * feq_common;

         const real_t t2x2 = t2x2_0 * rho;
         const real_t fac2 = t2x2 * inv2csq2;

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW;
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW;

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW;
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW;

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW;
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW;

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW;
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW;

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS;
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS;

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS;
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS;

         const real_t t1x2 = t1x2_0 * rho;
         const real_t fac1 = t1x2 * inv2csq2;

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S;
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S;

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W;
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W;

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B;
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B;

         const real_t t3x2 = t3x2_0 * rho;
         const real_t fac3 = t3x2 * inv2csq2;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         src->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW;
         src->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW;

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         src->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE;
         src->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE;

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         src->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW;
         src->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW;

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         src->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE;
         src->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE;
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6 )
{
   // relaxation parameter variables
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2 = real_t(2.0) / real_t(27.0) * real_t(2.0);  // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(54.0) * real_t(2.0);  // 1/54 * 2 else
   const real_t t3x2 = real_t(1.0) / real_t(216.0) * real_t(2.0); // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;
   const real_t fac3     = t3x2 * inv2csq2;

   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         dst->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common;// no force term

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         dst->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW + three_w2 * (  force[0] + force[1] );
         dst->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW - three_w2 * (  force[0] + force[1] );

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         dst->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW + three_w2 * (  force[0] - force[1] );
         dst->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW - three_w2 * (  force[0] - force[1] );

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         dst->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW + three_w2 * (  force[0] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW - three_w2 * (  force[0] + force[2] );

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         dst->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW + three_w2 * (  force[0] - force[2] );
         dst->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW - three_w2 * (  force[0] - force[2] );

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         dst->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS + three_w2 * (  force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS - three_w2 * (  force[1] + force[2] );

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         dst->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS + three_w2 * (  force[1] - force[2] );
         dst->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS - three_w2 * (  force[1] - force[2] );

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         dst->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S + three_w1 * force[1];
         dst->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S - three_w1 * force[1];

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         dst->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W + three_w1 * force[0];
         dst->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W - three_w1 * force[0];

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         dst->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B + three_w1 * force[2];
         dst->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B - three_w1 * force[2];

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         dst->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW + three_w3 * (  force[0] + force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW - three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         dst->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE + three_w3 * (  -force[0] + force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         dst->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW + three_w3 * (  force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         dst->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE + three_w3 * (  -force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6 )
{
   // relaxation parameter variables
   const real_t lambda_e = src->latticeModel().collisionModel().lambda_e();
   const real_t lambda_d = src->latticeModel().collisionModel().lambda_d();

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   // common prefactors for calculating the equilibrium parts
   const real_t t0   = real_t(8.0) / real_t(27.0);                // 8/27    for C
   const real_t t1x2 = real_t(2.0) / real_t(27.0) * real_t(2.0);  // 2/27 * 2 for N, S, W, E, T, B
   const real_t t2x2 = real_t(1.0) / real_t(54.0) * real_t(2.0);  // 1/54 * 2 else
   const real_t t3x2 = real_t(1.0) / real_t(216.0) * real_t(2.0); // 1/216    for TNE,BSW,TNW,BSE,TSE,BNW,TSW,BNE

   const real_t inv2csq2 = real_t(1.0) / ( real_t(2.0) * ( real_t(1.0) / real_t(3.0) ) * ( real_t(1.0) / real_t(3.0) ) ); //speed of sound related factor for equilibrium distribution function
   const real_t fac1     = t1x2 * inv2csq2;
   const real_t fac2     = t2x2 * inv2csq2;
   const real_t fac3     = t3x2 * inv2csq2;

   // relaxation parameter variables
   const real_t lambda_e_scaled = real_t(0.5) * lambda_e; // 0.5 times the usual value ...
   const real_t lambda_d_scaled = real_t(0.5) * lambda_d; // ... due to the way of calculations

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

         src->get( x, y, z, Stencil_T::idx[C] ) = vC * (real_t(1.0) - lambda_e) + lambda_e * t0 * feq_common; // no force term

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         const real_t velXPY = velX + velY;
         const real_t  sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common );
         const real_t asym_NE_SW = lambda_d_scaled * ( vNE - vSW - real_t(3.0) * t2x2 * velXPY );
         src->get( x, y, z, Stencil_T::idx[NE] ) = vNE - sym_NE_SW - asym_NE_SW + three_w2 * (  force[0] + force[1] );
         src->get( x, y, z, Stencil_T::idx[SW] ) = vSW - sym_NE_SW + asym_NE_SW - three_w2 * (  force[0] + force[1] );

         const real_t velXMY = velX - velY;
         const real_t  sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common );
         const real_t asym_SE_NW = lambda_d_scaled * ( vSE - vNW - real_t(3.0) * t2x2 * velXMY );
         src->get( x, y, z, Stencil_T::idx[SE] ) = vSE - sym_SE_NW - asym_SE_NW + three_w2 * (  force[0] - force[1] );
         src->get( x, y, z, Stencil_T::idx[NW] ) = vNW - sym_SE_NW + asym_SE_NW - three_w2 * (  force[0] - force[1] );

         const real_t velXPZ = velX + velZ;
         const real_t  sym_TE_BW = lambda_e_scaled * ( vTE + vBW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
         const real_t asym_TE_BW = lambda_d_scaled * ( vTE - vBW - real_t(3.0) * t2x2 * velXPZ );
         src->get( x, y, z, Stencil_T::idx[TE] ) = vTE - sym_TE_BW - asym_TE_BW + three_w2 * (  force[0] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BW] ) = vBW - sym_TE_BW + asym_TE_BW - three_w2 * (  force[0] + force[2] );

         const real_t velXMZ = velX - velZ;
         const real_t  sym_BE_TW = lambda_e_scaled * ( vBE + vTW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
         const real_t asym_BE_TW = lambda_d_scaled * ( vBE - vTW - real_t(3.0) * t2x2 * velXMZ );
         src->get( x, y, z, Stencil_T::idx[BE] ) = vBE - sym_BE_TW - asym_BE_TW + three_w2 * (  force[0] - force[2] );
         src->get( x, y, z, Stencil_T::idx[TW] ) = vTW - sym_BE_TW + asym_BE_TW - three_w2 * (  force[0] - force[2] );

         const real_t velYPZ = velY + velZ;
         const real_t  sym_TN_BS = lambda_e_scaled * ( vTN + vBS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
         const real_t asym_TN_BS = lambda_d_scaled * ( vTN - vBS - real_t(3.0) * t2x2 * velYPZ );
         src->get( x, y, z, Stencil_T::idx[TN] ) = vTN - sym_TN_BS - asym_TN_BS + three_w2 * (  force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BS] ) = vBS - sym_TN_BS + asym_TN_BS - three_w2 * (  force[1] + force[2] );

         const real_t velYMZ = velY - velZ;
         const real_t  sym_BN_TS = lambda_e_scaled * ( vBN + vTS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
         const real_t asym_BN_TS = lambda_d_scaled * ( vBN - vTS - real_t(3.0) * t2x2 * velYMZ );
         src->get( x, y, z, Stencil_T::idx[BN] ) = vBN - sym_BN_TS - asym_BN_TS + three_w2 * (  force[1] - force[2] );
         src->get( x, y, z, Stencil_T::idx[TS] ) = vTS - sym_BN_TS + asym_BN_TS - three_w2 * (  force[1] - force[2] );

         const real_t  sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common );
         const real_t asym_N_S = lambda_d_scaled * ( vN - vS - real_t(3.0) * t1x2 * velY );
         src->get( x, y, z, Stencil_T::idx[N] ) = vN - sym_N_S - asym_N_S + three_w1 * force[1];
         src->get( x, y, z, Stencil_T::idx[S] ) = vS - sym_N_S + asym_N_S - three_w1 * force[1];

         const real_t  sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common );
         const real_t asym_E_W = lambda_d_scaled * ( vE - vW - real_t(3.0) * t1x2 * velX );
         src->get( x, y, z, Stencil_T::idx[E] ) = vE - sym_E_W - asym_E_W + three_w1 * force[0];
         src->get( x, y, z, Stencil_T::idx[W] ) = vW - sym_E_W + asym_E_W - three_w1 * force[0];

         const real_t  sym_T_B = lambda_e_scaled * ( vT + vB  - fac1 * velZ * velZ - t1x2 * feq_common );
         const real_t asym_T_B = lambda_d_scaled * ( vT - vB - real_t(3.0) * t1x2 * velZ );
         src->get( x, y, z, Stencil_T::idx[T] ) = vT - sym_T_B - asym_T_B + three_w1 * force[2];
         src->get( x, y, z, Stencil_T::idx[B] ) = vB - sym_T_B + asym_T_B - three_w1 * force[2];

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t  sym_TNE_BSW = lambda_e_scaled * ( vTNE + vBSW - fac3 * vel_TNE_BSW * vel_TNE_BSW - t3x2 * feq_common );
         const real_t asym_TNE_BSW = lambda_d_scaled * ( vTNE - vBSW - real_t(3.0) * t3x2 * vel_TNE_BSW );
         src->get( x, y, z, Stencil_T::idx[TNE] ) = vTNE - sym_TNE_BSW - asym_TNE_BSW + three_w3 * (  force[0] + force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BSW] ) = vBSW - sym_TNE_BSW + asym_TNE_BSW - three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t  sym_TNW_BSE = lambda_e_scaled * ( vTNW + vBSE - fac3 * vel_TNW_BSE * vel_TNW_BSE - t3x2 * feq_common );
         const real_t asym_TNW_BSE = lambda_d_scaled * ( vTNW - vBSE - real_t(3.0) * t3x2 * vel_TNW_BSE );
         src->get( x, y, z, Stencil_T::idx[TNW] ) = vTNW - sym_TNW_BSE - asym_TNW_BSE + three_w3 * (  -force[0] + force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BSE] ) = vBSE - sym_TNW_BSE + asym_TNW_BSE - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t  sym_TSE_BNW = lambda_e_scaled * ( vTSE + vBNW - fac3 * vel_TSE_BNW * vel_TSE_BNW - t3x2 * feq_common );
         const real_t asym_TSE_BNW = lambda_d_scaled * ( vTSE - vBNW - real_t(3.0) * t3x2 * vel_TSE_BNW );
         src->get( x, y, z, Stencil_T::idx[TSE] ) = vTSE - sym_TSE_BNW - asym_TSE_BNW + three_w3 * (  force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNW] ) = vBNW - sym_TSE_BNW + asym_TSE_BNW - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t  sym_TSW_BNE = lambda_e_scaled * ( vTSW + vBNE - fac3 * vel_TSW_BNE * vel_TSW_BNE - t3x2 * feq_common );
         const real_t asym_TSW_BNE = lambda_d_scaled * ( vTSW - vBNE - real_t(3.0) * t3x2 * vel_TSW_BNE );
         src->get( x, y, z, Stencil_T::idx[TSW] ) = vTSW - sym_TSW_BNE - asym_TSW_BNE + three_w3 * (  -force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNE] ) = vBNE - sym_TSW_BNE + asym_TSW_BNE - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6' at the end of this file!






////////////////////////////////
// GENERIC TRT SPECIALIZATION //
////////////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT \
   std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value && \
   ! ( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7 )

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT )
{
   const real_t lambda_e =  lm.collisionModel().lambda_e();
   const real_t lambda_d =  lm.collisionModel().lambda_d();

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif

   real_t pdfs[ Stencil_T::Size ];

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP( src, numberOfGhostLayersToInclude, omp for schedule(static),

      if( this->filter(x,y,z) )
      {
         // stream pull
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const auto pdf = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
            dst->get( x, y, z, d.toIdx() ) = pdf;
            pdfs[ d.toIdx() ]              = pdf;
         }

         Vector3<real_t> velocity;
         real_t rho = this->densityVelocityIn( velocity, dst, x, y, z );

         if (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Guo_tag >::value)
            velocity -= real_c(0.5) * lm.forceModel().force(x,y,z);

         this->densityVelocityOut( x, y, z, lm, velocity, rho );

         const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho, lambda_e, lambda_e, lambda_d );

         // collide
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[ d.toIdx() ],
                                                                                           real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), lambda_e, lambda_e, lambda_d );

            const real_t fsym  = EquilibriumDistribution< LatticeModel_T >::getSymmetricPart ( *d, velocity, rho );
            const real_t fasym = EquilibriumDistribution< LatticeModel_T >::getAsymmetricPart( *d, velocity, rho );

            const real_t f     = pdfs[ d.toIdx() ];
            const real_t finv  = pdfs[ d.toInvIdx() ];

            dst->get( x, y, z, d.toIdx() ) = f - lambda_e * ( real_t( 0.5 ) * ( f + finv ) - fsym )
                                               - lambda_d * ( real_t( 0.5 ) * ( f - finv ) - fasym ) + forceTerm;
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP

#ifdef _OPENMP
   }
#endif
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT )
{
   const real_t lambda_e =  lm.collisionModel().lambda_e();
   const real_t lambda_d =  lm.collisionModel().lambda_d();

#ifdef _OPENMP
   #pragma omp parallel
   {
#endif

   real_t pdfs[ Stencil_T::Size ];

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP( src, numberOfGhostLayersToInclude, omp for schedule(static),

      if( this->filter(x,y,z) )
      {
         Vector3<real_t> velocity;
         real_t rho = this->densityVelocityIn( velocity, src, x, y, z );

         if (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Guo_tag >::value)
            velocity -= real_c(0.5) * lm.forceModel().force(x,y,z);

         this->densityVelocityOut( x, y, z, lm, velocity, rho );

         const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho, lambda_e, lambda_e, lambda_d );

         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
            pdfs[ d.toIdx() ] = src->get( x, y, z, d.toIdx() );

         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[ d.toIdx() ],
                                                                                           real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), lambda_e, lambda_e, lambda_d );

            const real_t fsym  = EquilibriumDistribution< LatticeModel_T >::getSymmetricPart ( *d, velocity, rho );
            const real_t fasym = EquilibriumDistribution< LatticeModel_T >::getAsymmetricPart( *d, velocity, rho );

            const real_t f     = pdfs[ d.toIdx() ];
            const real_t finv  = pdfs[ d.toInvIdx() ];

            src->get( x, y, z, d.toIdx() ) = f - lambda_e * ( real_t( 0.5 ) * ( f + finv ) - fsym )
                                               - lambda_d * ( real_t( 0.5 ) * ( f - finv ) - fasym ) + forceTerm;
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP

#ifdef _OPENMP
   }
#endif
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_1
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_2
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_3
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_4
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_5
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_6
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_TRT_7



} // namespace lbm
} // namespace walberla
