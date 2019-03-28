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

// For a generic implementation of the SRT sweep jump to the end of this file!



//////////////////////////
// D2Q9 SPECIALIZATIONS //
//////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9 \
    (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
    std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value && \
    LatticeModel_T::CollisionModel::constant && \
    ! LatticeModel_T::compressible && \
    LatticeModel_T::equilibriumAccuracyOrder == 2 && \
    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
    std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(4) / real_t( 9) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t( 9) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, real_t(0) ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(4) / real_t( 9) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t( 9) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, real_t(0) ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9' at the end of this file!






///////////////////////////
// D3Q19 SPECIALIZATIONS //
///////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );
         
         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ );
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ );
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm;

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX );
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX );
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY );
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY );
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ );
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ );

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm;

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX );
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX );
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY );
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY );
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ );
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ );

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2' at the end of this file!





             
#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm; // no force term

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;
         
         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm; // no force term

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;
         
         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm; // no force term

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;
         
         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(1) / real_t( 3) ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(1) / real_t(18) ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1) / real_t(36) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(1) / real_t(6) );
   const real_t three_w2( real_t(1) / real_t(12) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP()
         
         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm; // no force term

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;
         
         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4' at the end of this file!






///////////////////////////
// D3Q27 SPECIALIZATIONS //
///////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ );
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ );

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         dst->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3 * ( vel_trm_TNE_BSW + vel_TNE_BSW );
         dst->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3 * ( vel_trm_TNE_BSW - vel_TNE_BSW );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         dst->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3 * ( vel_trm_TNW_BSE + vel_TNW_BSE );
         dst->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3 * ( vel_trm_TNW_BSE - vel_TNW_BSE );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         dst->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3 * ( vel_trm_TSE_BNW + vel_TSE_BNW );
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3 * ( vel_trm_TSE_BNW - vel_TSE_BNW );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         dst->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3 * ( vel_trm_TSW_BNE + vel_TSW_BNE );
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3 * ( vel_trm_TSW_BNE - vel_TSW_BNE );

      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX );
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX );
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY );
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY );
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ );
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ );

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ );

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         src->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3 * ( vel_trm_TNE_BSW + vel_TNE_BSW );
         src->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3 * ( vel_trm_TNE_BSW - vel_TNE_BSW );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         src->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3 * ( vel_trm_TNW_BSE + vel_TNW_BSE );
         src->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3 * ( vel_trm_TNW_BSE - vel_TNW_BSE );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         src->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3 * ( vel_trm_TSE_BNW + vel_TSE_BNW );
         src->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3 * ( vel_trm_TSE_BNW - vel_TSE_BNW );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         src->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3 * ( vel_trm_TSW_BNE + vel_TSW_BNE );
         src->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3 * ( vel_trm_TSW_BNE - vel_TSW_BNE );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm;

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX );
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX );
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY );
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY );
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ );
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ );

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );

         const real_t omega_w3_rho = omega_w3 * rho;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         dst->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3_rho * ( vel_trm_TNE_BSW + vel_TNE_BSW );
         dst->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3_rho * ( vel_trm_TNE_BSW - vel_TNE_BSW );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         dst->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3_rho * ( vel_trm_TNW_BSE + vel_TNW_BSE );
         dst->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3_rho * ( vel_trm_TNW_BSE - vel_TNW_BSE );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         dst->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3_rho * ( vel_trm_TSE_BNW + vel_TSE_BNW );
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3_rho * ( vel_trm_TSE_BNW - vel_TSE_BNW );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         dst->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3_rho * ( vel_trm_TSW_BNE + vel_TSW_BNE );
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3_rho * ( vel_trm_TSW_BNE - vel_TSW_BNE );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm;

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX );
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX );
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY );
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY );
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ );
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ );

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );

         const real_t omega_w3_rho = omega_w3 * rho;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         src->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3_rho * ( vel_trm_TNE_BSW + vel_TNE_BSW );
         src->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3_rho * ( vel_trm_TNE_BSW - vel_TNE_BSW );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         src->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3_rho * ( vel_trm_TNW_BSE + vel_TNW_BSE );
         src->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3_rho * ( vel_trm_TNW_BSE - vel_TNW_BSE );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         src->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3_rho * ( vel_trm_TSE_BNW + vel_TSE_BNW );
         src->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3_rho * ( vel_trm_TSE_BNW - vel_TSE_BNW );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         src->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3_rho * ( vel_trm_TSW_BNE + vel_TSW_BNE );
         src->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3_rho * ( vel_trm_TSW_BNE - vel_TSW_BNE );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm; // no force term

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         dst->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3 * ( vel_trm_TNE_BSW + vel_TNE_BSW )+ three_w3 * (  force[0] + force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3 * ( vel_trm_TNE_BSW - vel_TNE_BSW )- three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         dst->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3 * ( vel_trm_TNW_BSE + vel_TNW_BSE ) + three_w3 * (  -force[0] + force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3 * ( vel_trm_TNW_BSE - vel_TNW_BSE ) - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         dst->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3 * ( vel_trm_TSE_BNW + vel_TSE_BNW ) + three_w3 * (  force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3 * ( vel_trm_TSE_BNW - vel_TSE_BNW ) - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         dst->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3 * ( vel_trm_TSW_BNE + vel_TSW_BNE ) + three_w3 * (  -force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3 * ( vel_trm_TSW_BNE - vel_TSW_BNE ) - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho + real_t(1) );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third * rho - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * dir_indep_trm; // no force term

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1 * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1 * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2 * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2 * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2 * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2 * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2 * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2 * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2 * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2 * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2 * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         src->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3 * ( vel_trm_TNE_BSW + vel_TNE_BSW ) + three_w3 * (  force[0] + force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3 * ( vel_trm_TNE_BSW - vel_TNE_BSW ) - three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         src->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3 * ( vel_trm_TNW_BSE + vel_TNW_BSE ) + three_w3 * (  -force[0] + force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3 * ( vel_trm_TNW_BSE - vel_TNW_BSE ) - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         src->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3 * ( vel_trm_TSE_BNW + vel_TSE_BNW ) + three_w3 * (  force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3 * ( vel_trm_TSE_BNW - vel_TSE_BNW ) - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         src->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3 * ( vel_trm_TSW_BNE + vel_TSW_BNE ) + three_w3 * (  -force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3 * ( vel_trm_TSW_BNE - vel_TSW_BNE ) - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7' at the end of this file!






#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::CollisionModel::constant && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Simple_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         dst->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm; // no force term

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         dst->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         dst->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         dst->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         dst->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         dst->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         dst->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         dst->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         dst->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         dst->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         dst->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         dst->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         dst->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         dst->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );

         const real_t omega_w3_rho = omega_w3 * rho;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         dst->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3_rho * ( vel_trm_TNE_BSW + vel_TNE_BSW ) + three_w3 * (  force[0] + force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3_rho * ( vel_trm_TNE_BSW - vel_TNE_BSW ) - three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         dst->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3_rho * ( vel_trm_TNW_BSE + vel_TNW_BSE ) + three_w3 * (  -force[0] + force[1] + force[2] );
         dst->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3_rho * ( vel_trm_TNW_BSE - vel_TNW_BSE ) - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         dst->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3_rho * ( vel_trm_TSE_BNW + vel_TSE_BNW ) + three_w3 * (  force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3_rho * ( vel_trm_TSE_BNW - vel_TSE_BNW ) - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         dst->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3_rho * ( vel_trm_TSW_BNE + vel_TSW_BNE ) + three_w3 * (  -force[0] - force[1] + force[2] );
         dst->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3_rho * ( vel_trm_TSW_BNE - vel_TSW_BNE ) - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8 )
{
   const real_t omega = src->latticeModel().collisionModel().omega();

   const real_t omega_trm( real_t(1) - omega );
   const real_t  omega_w0( real_t(3) * ( real_t(8.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w1( real_t(3) * ( real_t(2.0) / real_t(27.0)  ) * omega );
   const real_t  omega_w2( real_t(3) * ( real_t(1.0) / real_t(54.0)  ) * omega );
   const real_t  omega_w3( real_t(3) * ( real_t(1.0) / real_t(216.0) ) * omega );
   const real_t one_third( real_t(1) / real_t(3) );

   const real_t three_w1( real_t(2) / real_t(9) );
   const real_t three_w2( real_t(1) / real_t(18) );
   const real_t three_w3( real_t(1) / real_t(72) );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), rho );

         const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;

         const real_t dir_indep_trm = one_third - real_t(0.5) * ( velXX + velYY + velZZ );

         src->get(x,y,z,Stencil_T::idx[C]) = omega_trm * vC + omega_w0 * rho * dir_indep_trm; // no force term

         const real_t omega_w1_rho = omega_w1 * rho;

         const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
         const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
         const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

         const Vector3< real_t > & force = src->latticeModel().forceModel().force(x,y,z);

         src->get(x,y,z,Stencil_T::idx[E]) = omega_trm * vE + omega_w1_rho * ( vel_trm_E_W + velX ) + three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[W]) = omega_trm * vW + omega_w1_rho * ( vel_trm_E_W - velX ) - three_w1 * force[0];
         src->get(x,y,z,Stencil_T::idx[N]) = omega_trm * vN + omega_w1_rho * ( vel_trm_N_S + velY ) + three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[S]) = omega_trm * vS + omega_w1_rho * ( vel_trm_N_S - velY ) - three_w1 * force[1];
         src->get(x,y,z,Stencil_T::idx[T]) = omega_trm * vT + omega_w1_rho * ( vel_trm_T_B + velZ ) + three_w1 * force[2];
         src->get(x,y,z,Stencil_T::idx[B]) = omega_trm * vB + omega_w1_rho * ( vel_trm_T_B - velZ ) - three_w1 * force[2];

         const real_t omega_w2_rho = omega_w2 * rho;

         const real_t velXmY = velX - velY;
         const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

         src->get(x,y,z,Stencil_T::idx[NW]) = omega_trm * vNW + omega_w2_rho * ( vel_trm_NW_SE - velXmY ) + three_w2 * (  force[1] - force[0] );
         src->get(x,y,z,Stencil_T::idx[SE]) = omega_trm * vSE + omega_w2_rho * ( vel_trm_NW_SE + velXmY ) + three_w2 * (  force[0] - force[1] );

         const real_t velXpY = velX + velY;
         const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

         src->get(x,y,z,Stencil_T::idx[NE]) = omega_trm * vNE + omega_w2_rho * ( vel_trm_NE_SW + velXpY ) + three_w2 * (  force[0] + force[1] );
         src->get(x,y,z,Stencil_T::idx[SW]) = omega_trm * vSW + omega_w2_rho * ( vel_trm_NE_SW - velXpY ) + three_w2 * ( -force[0] - force[1] );

         const real_t velXmZ = velX - velZ;
         const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

         src->get(x,y,z,Stencil_T::idx[TW]) = omega_trm * vTW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ ) + three_w2 * (  force[2] - force[0] );
         src->get(x,y,z,Stencil_T::idx[BE]) = omega_trm * vBE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ ) + three_w2 * (  force[0] - force[2] );

         const real_t velXpZ = velX + velZ;
         const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

         src->get(x,y,z,Stencil_T::idx[TE]) = omega_trm * vTE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ ) + three_w2 * (  force[0] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BW]) = omega_trm * vBW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ ) + three_w2 * ( -force[0] - force[2] );

         const real_t velYmZ = velY - velZ;
         const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

         src->get(x,y,z,Stencil_T::idx[TS]) = omega_trm * vTS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ ) + three_w2 * (  force[2] - force[1] );
         src->get(x,y,z,Stencil_T::idx[BN]) = omega_trm * vBN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ ) + three_w2 * (  force[1] - force[2] );

         const real_t velYpZ = velY + velZ;
         const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

         src->get(x,y,z,Stencil_T::idx[TN]) = omega_trm * vTN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ ) + three_w2 * (  force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BS]) = omega_trm * vBS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ ) + three_w2 * ( -force[1] - force[2] );

         const real_t omega_w3_rho = omega_w3 * rho;

         const real_t vel_TNE_BSW = velX + velY + velZ;
         const real_t vel_trm_TNE_BSW = dir_indep_trm + real_t(1.5) * vel_TNE_BSW * vel_TNE_BSW;

         src->get(x,y,z,Stencil_T::idx[TNE]) = omega_trm * vTNE + omega_w3_rho * ( vel_trm_TNE_BSW + vel_TNE_BSW ) + three_w3 * (  force[0] + force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BSW]) = omega_trm * vBSW + omega_w3_rho * ( vel_trm_TNE_BSW - vel_TNE_BSW ) - three_w3 * (  force[0] + force[1] + force[2] );

         const real_t vel_TNW_BSE = -velX + velY + velZ;
         const real_t vel_trm_TNW_BSE = dir_indep_trm + real_t(1.5) * vel_TNW_BSE * vel_TNW_BSE;

         src->get(x,y,z,Stencil_T::idx[TNW]) = omega_trm * vTNW + omega_w3_rho * ( vel_trm_TNW_BSE + vel_TNW_BSE ) + three_w3 * (  -force[0] + force[1] + force[2] );
         src->get(x,y,z,Stencil_T::idx[BSE]) = omega_trm * vBSE + omega_w3_rho * ( vel_trm_TNW_BSE - vel_TNW_BSE ) - three_w3 * (  -force[0] + force[1] + force[2] );

         const real_t vel_TSE_BNW = velX - velY + velZ;
         const real_t vel_trm_TSE_BNW = dir_indep_trm + real_t(1.5) * vel_TSE_BNW * vel_TSE_BNW;

         src->get( x, y, z, Stencil_T::idx[TSE] ) = omega_trm * vTSE + omega_w3_rho * ( vel_trm_TSE_BNW + vel_TSE_BNW ) + three_w3 * (  force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNW] ) = omega_trm * vBNW + omega_w3_rho * ( vel_trm_TSE_BNW - vel_TSE_BNW ) - three_w3 * (  force[0] - force[1] + force[2] );

         const real_t vel_TSW_BNE = - velX - velY + velZ;
         const real_t vel_trm_TSW_BNE = dir_indep_trm + real_t(1.5) * vel_TSW_BNE * vel_TSW_BNE;

         src->get( x, y, z, Stencil_T::idx[TSW] ) = omega_trm * vTSW + omega_w3_rho * ( vel_trm_TSW_BNE + vel_TSW_BNE ) + three_w3 * (  -force[0] - force[1] + force[2] );
         src->get( x, y, z, Stencil_T::idx[BNE] ) = omega_trm * vBNE + omega_w3_rho * ( vel_trm_TSW_BNE - vel_TSW_BNE ) - three_w3 * (  -force[0] - force[1] + force[2] );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// Do not forget to exclude 'WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8' from the generic implementation at the end of the file!
// Do not forget the '#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8' at the end of this file!






////////////////////////////////
// GENERIC SRT SPECIALIZATION //
////////////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT \
   std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value && \
   ! ( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8 || \
       WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9 )

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT )
{
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         // stream pull
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
            dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

         Vector3<real_t> velocity;
         real_t rho = this->densityVelocityIn( velocity, dst, x, y, z );

         this->densityVelocityOut( x, y, z, lm, velocity, rho );

         const real_t omega = lm.collisionModel().omega( x, y, z, velocity, rho );

         const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho, omega, omega );

         // collide
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[ d.toIdx() ],
                                                                                           real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), omega, omega );

            dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega ) * dst->get( x, y, z, d.toIdx() ) +
                                                             omega   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                                             forceTerm;
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT )
{
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         Vector3<real_t> velocity;
         real_t rho = this->densityVelocityIn( velocity, src, x, y, z );

         this->densityVelocityOut( x, y, z, lm, velocity, rho );

         const real_t omega = lm.collisionModel().omega( x, y, z, velocity, rho );

         const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho, omega, omega );

         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[ d.toIdx() ],
                                                                                           real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), omega, omega );

            src->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega ) * src->get( x, y, z, d.toIdx() ) +
                                                             omega   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                                             forceTerm;
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_1
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_2
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_3
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_4
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_5
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_6
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_7
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_8
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_SRT_9



} // namespace lbm
} // namespace walberla
