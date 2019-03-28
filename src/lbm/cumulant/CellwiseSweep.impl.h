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
//! \author Sagar Dolas <sagar.dolas@fau.de>
//! \author Sunil Konatham <sunil.Konatham@fau.de>

//======================================================================================================================


#include "core/math/FPClassify.h"

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


namespace walberla{
namespace lbm{

///////////////////////////
// D3Q27 SPECIALIZATIONS //
///////////////////////////

// I will need to some changes to this also 
#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_CUMULANT_1 \
   (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::Cumulant_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value && \
   LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value)
   
WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_CUMULANT_1 ) // changed to SRT_1 for cumulant only in this file 

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_CUMULANT_1 )
{
  //WALBERLA_LOG_DEVEL("Cumulant StreamCollide");
  
   // retrieving all omegas
   const real_t omega1 = src->latticeModel().collisionModel().omega1(); 
   const real_t omega2 = src->latticeModel().collisionModel().omega2();
   const real_t omega3 = src->latticeModel().collisionModel().omega3();
   const real_t omega4 = src->latticeModel().collisionModel().omega4();
   const real_t omega5 = src->latticeModel().collisionModel().omega5();
   const real_t omega6 = src->latticeModel().collisionModel().omega6();
   const real_t omega7 = src->latticeModel().collisionModel().omega7();
   const real_t omega8 = src->latticeModel().collisionModel().omega8();
   const real_t omega9 = src->latticeModel().collisionModel().omega9();
   const real_t omega10 = src->latticeModel().collisionModel().omega10();
   
   const real_t omega_trm1( real_t(1.0) - omega1 );
   const real_t omega_trm2( real_t(1.0) - omega2 );
   const real_t omega_trm3( real_t(1.0) - omega3 );
   const real_t omega_trm4( real_t(1.0) - omega4 );
   const real_t omega_trm5( real_t(1.0) - omega5 );
   const real_t omega_trm6( real_t(1.0) - omega6 );
   const real_t omega_trm7( real_t(1.0) - omega7 );
   const real_t omega_trm8( real_t(1.0) - omega8 );
   const real_t omega_trm9( real_t(1.0) - omega9 );
   const real_t omega_trm10( real_t(1.0) - omega10 );
  
   
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) ) 
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP() 
	 
	 const real_t updated_rho = rho;

         this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), updated_rho );  // dont really understand what is been done here 
	 
         // defining the sqaure of velocities 
	 const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;
	 
	 // defining velocity terms for central moment space 
	 
	 // velocity in z direction 
	 const real_t velz_term = real_t(-1.0) * velZ ;
	 const real_t sqr_velz_term = velz_term * velz_term ;
	 const real_t oneminus_velz = real_t(1.0) - velZ ;
	 const real_t sqr_oneminus_velz = oneminus_velz * oneminus_velz ;
	 const real_t negoneminus_velz = real_t(-1.0) - velZ ;
	 const real_t sqr_negoneminus_velz = negoneminus_velz * negoneminus_velz ; 
	 
	 // Declaring the constants and initialising them with distribution function to transform  for the moment space 
	 const real_t k_00_0 = vBSW                        + vSW                    + vTSW ;
	 const real_t k_00_1 = vBSW * negoneminus_velz     + vSW * velz_term        + vTSW * oneminus_velz ;
	 const real_t k_00_2 = vBSW * sqr_negoneminus_velz + vSW * sqr_velz_term    + vTSW * sqr_oneminus_velz ;
	 const real_t k_01_0 = vBW                         + vW                     + vTW  ; 
	 const real_t k_01_1 = vBW * negoneminus_velz      + vW  * velz_term        + vTW * oneminus_velz ;
	 const real_t k_01_2 = vBW * sqr_negoneminus_velz  + vW  * sqr_velz_term    + vTW * sqr_oneminus_velz ;
	 const real_t k_02_0 = vBNW                        + vNW                    + vTNW ;
	 const real_t k_02_1 = vBNW * negoneminus_velz     + vNW * velz_term        + vTNW * oneminus_velz ;
	 const real_t k_02_2 = vBNW * sqr_negoneminus_velz + vNW * sqr_velz_term    + vTNW * sqr_oneminus_velz;
	 const real_t k_10_0 = vBS                         + vS                     + vTS;
	 const real_t k_10_1 = vBS * negoneminus_velz      + vS * velz_term         + vTS * oneminus_velz ;
	 const real_t k_10_2 = vBS * sqr_negoneminus_velz  + vS * sqr_velz_term     + vTS * sqr_oneminus_velz ;
	 const real_t k_11_0 = vB                          + vC                     + vT ;
	 const real_t k_11_1 = vB * negoneminus_velz       + vC * velz_term         + vT * oneminus_velz ;
	 const real_t k_11_2 = vB * sqr_negoneminus_velz   + vC * sqr_velz_term     + vT * sqr_oneminus_velz ;
	 const real_t k_12_0 = vBN                         + vN                     + vTN ;
	 const real_t k_12_1 = vBN * negoneminus_velz      + vN * velz_term         + vTN * oneminus_velz ;
	 const real_t k_12_2 = vBN * sqr_negoneminus_velz  + vN * sqr_velz_term     + vTN * sqr_oneminus_velz ;
	 const real_t k_20_0 = vBSE                        + vSE                    + vTSE ;
	 const real_t k_20_1 = vBSE * negoneminus_velz     + vSE * velz_term        + vTSE * oneminus_velz ;
	 const real_t k_20_2 = vBSE * sqr_negoneminus_velz + vSE * sqr_velz_term    + vTSE * sqr_oneminus_velz ;
	 const real_t k_21_0 = vBE                         + vE                     + vTE ;
	 const real_t k_21_1 = vBE * negoneminus_velz      + vE * velz_term         + vTE * oneminus_velz ;
	 const real_t k_21_2 = vBE * sqr_negoneminus_velz  + vE * sqr_velz_term     + vTE * sqr_oneminus_velz ;
	 const real_t k_22_0 = vBNE                        + vNE                    + vTNE ;
	 const real_t k_22_1 = vBNE * negoneminus_velz     + vNE * velz_term        + vTNE * oneminus_velz ;
	 const real_t k_22_2 = vBNE * sqr_negoneminus_velz + vNE * sqr_velz_term    + vTNE * sqr_oneminus_velz ;
	 
	 // velocity in y direction 
	 const real_t vely_term = real_t(-1.0) * velY ;
	 const real_t sqr_vely_term = vely_term * vely_term ;
	 const real_t oneminus_vely = real_t(1.0) - velY ;
	 const real_t sqr_oneminus_vely = oneminus_vely * oneminus_vely ;
	 const real_t negoneminus_vely = real_t(-1.0) - velY ;
	 const real_t sqr_negoneminus_vely = negoneminus_vely * negoneminus_vely ;
	 
	 const real_t k_0_00 = k_00_0				+ k_01_0			+ k_02_0;
	 const real_t k_0_10 = k_00_0 * negoneminus_vely	+ k_01_0 * vely_term		+ k_02_0 * oneminus_vely;
	 const real_t k_0_20 = k_00_0 * sqr_negoneminus_vely	+ k_01_0 * sqr_vely_term	+ k_02_0 * sqr_oneminus_vely; 
	 const real_t k_0_01 = k_00_1				+ k_01_1			+ k_02_1;
	 const real_t k_0_11 = k_00_1 * negoneminus_vely	+ k_01_1 * vely_term		+ k_02_1 *  oneminus_vely;
	 const real_t k_0_21 = k_00_1 * sqr_negoneminus_vely	+ k_01_1 * sqr_vely_term	+ k_02_1 * sqr_oneminus_vely;
	 const real_t k_0_02 = k_00_2				+ k_01_2			+ k_02_2;
	 const real_t k_0_12 = k_00_2 * negoneminus_vely	+ k_01_2 * vely_term		+ k_02_2 *  oneminus_vely;
	 const real_t k_0_22 = k_00_2 * sqr_negoneminus_vely	+ k_01_2 * sqr_vely_term	+ k_02_2 * sqr_oneminus_vely;
	 const real_t k_1_00 = k_10_0				+ k_11_0			+ k_12_0;
	 const real_t k_1_10 = k_10_0 * negoneminus_vely	+ k_11_0 * vely_term		+ k_12_0 *  oneminus_vely;
	 const real_t k_1_20 = k_10_0 * sqr_negoneminus_vely	+ k_11_0 * sqr_vely_term	+ k_12_0 * sqr_oneminus_vely;
	 const real_t k_1_01 = k_10_1				+ k_11_1			+ k_12_1;
	 const real_t k_1_11 = k_10_1 * negoneminus_vely	+ k_11_1 * vely_term		+ k_12_1 *  oneminus_vely;
	 const real_t k_1_21 = k_10_1 * sqr_negoneminus_vely	+ k_11_1 * sqr_vely_term	+ k_12_1 * sqr_oneminus_vely;
	 const real_t k_1_02 = k_10_2 				+ k_11_2			+ k_12_2 ;
	 const real_t k_1_12 = k_10_2 * negoneminus_vely	+ k_11_2 * vely_term		+ k_12_2 *  oneminus_vely;
	 const real_t k_1_22 = k_10_2 * sqr_negoneminus_vely	+ k_11_2 * sqr_vely_term	+ k_12_2 * sqr_oneminus_vely;
	 const real_t k_2_00 = k_20_0				+ k_21_0			+ k_22_0;
	 const real_t k_2_10 = k_20_0 * negoneminus_vely	+ k_21_0 * vely_term		+ k_22_0 *  oneminus_vely;
	 const real_t k_2_20 = k_20_0 * sqr_negoneminus_vely	+ k_21_0 * sqr_vely_term	+ k_22_0 * sqr_oneminus_vely;
	 const real_t k_2_01 = k_20_1				+ k_21_1			+ k_22_1;
	 const real_t k_2_11 = k_20_1 * negoneminus_vely	+ k_21_1 * vely_term		+ k_22_1 *  oneminus_vely;
	 const real_t k_2_21 = k_20_1 * sqr_negoneminus_vely	+ k_21_1 * sqr_vely_term	+ k_22_1 * sqr_oneminus_vely;
	 const real_t k_2_02 = k_20_2				+ k_21_2			+ k_22_2;
	 const real_t k_2_12 = k_20_2 * negoneminus_vely	+ k_21_2 * vely_term		+ k_22_2 *  oneminus_vely;
	 const real_t k_2_22 = k_20_2 * sqr_negoneminus_vely	+ k_21_2 * sqr_vely_term	+ k_22_2 * sqr_oneminus_vely;
	 
	 // Velocity in x direction 
	 const real_t velx_term = real_t(-1.0) * velX ;
	 const real_t sqr_velx_term = velx_term * velx_term ;
	 const real_t oneminus_velx = real_t(1.0) - velX ;
	 const real_t sqr_oneminus_velx = oneminus_velx * oneminus_velx ;
	 const real_t negoneminus_velx = real_t(-1.0) - velX ;
	 const real_t sqr_negoneminus_velx = negoneminus_velx * negoneminus_velx ;
	 
	 const real_t k_000 = k_0_00				+ k_1_00			+ k_2_00;
	 const real_t k_100 = k_0_00 * negoneminus_velx		+ k_1_00 * velx_term		+ k_2_00 * oneminus_velx;
	 const real_t k_200 = k_0_00 * sqr_negoneminus_velx	+ k_1_00 * sqr_velx_term	+ k_2_00 * sqr_oneminus_velx;
	 const real_t k_001 = k_0_01				+ k_1_01			+ k_2_01;
	 const real_t k_101 = k_0_01 * negoneminus_velx		+ k_1_01 * velx_term		+ k_2_01 * oneminus_velx;
	 const real_t k_201 = k_0_01 * sqr_negoneminus_velx	+ k_1_01 * sqr_velx_term	+ k_2_01 * sqr_oneminus_velx;
	 const real_t k_002 = k_0_02				+ k_1_02			+ k_2_02;
	 const real_t k_102 = k_0_02 * negoneminus_velx		+ k_1_02 * velx_term		+ k_2_02 * oneminus_velx;
	 const real_t k_202 = k_0_02 * sqr_negoneminus_velx	+ k_1_02 * sqr_velx_term	+ k_2_02 * sqr_oneminus_velx;
	 const real_t k_010 = k_0_10				+ k_1_10			+ k_2_10;
	 const real_t k_110 = k_0_10 * negoneminus_velx		+ k_1_10 * velx_term		+ k_2_10 * oneminus_velx;
	 const real_t k_210 = k_0_10 * sqr_negoneminus_velx	+ k_1_10 * sqr_velx_term	+ k_2_10 * sqr_oneminus_velx;
	 const real_t k_011 = k_0_11				+ k_1_11			+ k_2_11;
	 const real_t k_111 = k_0_11 * negoneminus_velx		+ k_1_11 * velx_term		+ k_2_11 * oneminus_velx;
	 const real_t k_211 = k_0_11 * sqr_negoneminus_velx	+ k_1_11 * sqr_velx_term	+ k_2_11 * sqr_oneminus_velx;
	 const real_t k_012 = k_0_12				+ k_1_12			+ k_2_12;
	 const real_t k_112 = k_0_12 * negoneminus_velx		+ k_1_12 * velx_term		+ k_2_12 * oneminus_velx;
	 const real_t k_212 = k_0_12 * sqr_negoneminus_velx	+ k_1_12 * sqr_velx_term	+ k_2_12 * sqr_oneminus_velx;
	 const real_t k_020 = k_0_20				+ k_1_20			+ k_2_20;
	 const real_t k_120 = k_0_20 * negoneminus_velx		+ k_1_20 * velx_term		+ k_2_20 * oneminus_velx;
	 const real_t k_220 = k_0_20 * sqr_negoneminus_velx	+ k_1_20 * sqr_velx_term	+ k_2_20 * sqr_oneminus_velx;
	 const real_t k_021 = k_0_21				+ k_1_21			+ k_2_21;
	 const real_t k_121 = k_0_21 * negoneminus_velx		+ k_1_21 * velx_term		+ k_2_21 * oneminus_velx;
	 const real_t k_221 = k_0_21 * sqr_negoneminus_velx	+ k_1_21 * sqr_velx_term	+ k_2_21 * sqr_oneminus_velx;
	 const real_t k_022 = k_0_22				+ k_1_22			+ k_2_22;
	 const real_t k_122 = k_0_22 * negoneminus_velx		+ k_1_22 * velx_term		+ k_2_22 * oneminus_velx;
	 const real_t k_222 = k_0_22 * sqr_negoneminus_velx	+ k_1_22 * sqr_velx_term	+ k_2_22 * sqr_oneminus_velx;
	 
	 // transformation to moment space done now , transform to cumulant space 
	 
	 // defining the constants for central moment space 
	 const real_t rho_inv = real_t(1.0) / updated_rho; 

	 // defining the sqaures 
	 const real_t sqr_k_110 = k_110 * k_110;
	 const real_t sqr_k_101 = k_101 * k_101;
	 const real_t sqr_k_011 = k_011 * k_011;
	 const real_t sqr_k_111 = k_111 * k_111;
	 
	 const real_t sqr_rho_inv = rho_inv * rho_inv ;
	 
	 // defining the cumulant space constants
	 const real_t C_000 = k_000 ;
	 const real_t C_001 = k_001 ;
	 const real_t C_002 = k_002 ;
	 const real_t C_010 = k_010 ;
	 const real_t C_011 = k_011 ;
	 const real_t C_012 = k_012 ;
	 const real_t C_020 = k_020 ;
	 const real_t C_021 = k_021 ;
	 const real_t C_022 = k_022 - (k_020 * k_002 + real_t(2.0) * sqr_k_011) * rho_inv ;
	 const real_t C_100 = k_100 ;
	 const real_t C_101 = k_101 ;
	 const real_t C_102 = k_102 ;
	 const real_t C_110 = k_110 ;
	 const real_t C_111 = k_111 ;
	 const real_t C_112 = k_112 - (k_002 * k_110 + real_t(2.0) * k_101 * k_011) * rho_inv ;
	 const real_t C_120 = k_120 ;
	 const real_t C_121 = k_121 - (k_020 * k_101 + real_t(2.0) * k_011 * k_110) * rho_inv ;
	 const real_t C_122 = k_122 - (k_002 * k_120 + k_020 * k_120 + real_t(4.0) * k_011 * k_111 + real_t(2.0) * (k_101 * k_021 + k_110 * k_012)) * rho_inv ;
	 const real_t C_200 = k_200 ;
	 const real_t C_201 = k_201 ;
	 const real_t C_202 = k_202 - (k_002 * k_200 + real_t(2.0) * sqr_k_101) * rho_inv ;
	 const real_t C_210 = k_210 ;
	 const real_t C_211 = k_211 - (k_200 * k_011 + real_t(2.0) * k_110 * k_101) * rho_inv ;
	 const real_t C_212 = k_212 - (k_002 * k_210 + k_200 * k_210 + real_t(4.0) * k_101 * k_111 + real_t(2.0) * (k_011 * k_201 + k_110 * k_102)) * rho_inv ;
	 const real_t C_220 = k_220 - (k_200 * k_020 + real_t(2.0) * sqr_k_110) * rho_inv ;
	 const real_t C_221 = k_221 - (k_020 * k_201 + k_200 * k_201 + real_t(4.0) * k_110 * k_111 + real_t(2.0) * (k_011 * k_210 + k_101 * k_120)) * rho_inv;
	 const real_t C_222 = k_222 - (real_t(4.0)* sqr_k_111 + k_200 * k_022 + k_020 * k_202 + k_002 * k_220  + real_t(4.0) * (k_011 * k_211 + k_101 * k_121 + k_110 * k_112) + real_t(2.0) * (k_120 * k_102 + k_210 * k_012 + k_201 * k_021)) * rho_inv 
 				    + (real_t(16.0) * k_110 * k_101 * k_011 + real_t(4.0) * (k_020 * sqr_k_101 + k_200 * sqr_k_011 + k_002 * sqr_k_110) + real_t(2.0) * k_200 *k_020 * k_002) * sqr_rho_inv;
	 
	 
	 // collision happens in cumulant space 
	 const real_t CS_000 = C_000 ;
	 const real_t CS_100 = C_100 ;
	 const real_t CS_010 = C_010 ;
	 const real_t CS_001 = C_001 ;
	 const real_t CS_110 = omega_trm1 * C_110 ;	
	 const real_t CS_101 = omega_trm1 * C_101 ;
	 const real_t CS_011 = omega_trm1 * C_011 ;
	 
	 const real_t Dxux = real_t(-0.5) * omega1 * rho_inv * (real_t(2.0) * C_200 - C_020 - C_002) - real_t(0.5) * omega2 * rho_inv * (C_200 + C_020 + C_002 - C_000);
	 const real_t Dyuy = Dxux + real_t(1.5) * omega1 * rho_inv * (C_200 - C_020);
	 const real_t Dzuz = Dxux + real_t(1.5) * omega1 * rho_inv * (C_200 - C_002);
	 
	 const real_t CS_200__m__CS020 = (omega_trm1) * (C_200 - C_020) - real_t(3.0) * updated_rho * (real_t(1.0) - real_t(0.5) * omega1) * (Dxux * velXX - Dyuy * velYY) ;
	 const real_t CS_200__m__CS002 = (omega_trm1) * (C_200 - C_002) - real_t(3.0) * updated_rho * (real_t(1.0) - real_t(0.5) * omega1) * (Dxux * velXX - Dzuz * velZZ) ;
	 const real_t CS_200__p__CS020__p__CS_002 = omega2 * C_000 + (omega_trm2) * (C_200 + C_020 + C_002) - real_t(3.0) * updated_rho *(real_t(1.0) - real_t(0.5) * omega2) * (Dxux * velXX + Dyuy * velYY + Dzuz * velZZ);
	 const real_t CS_200 = (CS_200__m__CS020 + CS_200__m__CS002 + CS_200__p__CS020__p__CS_002) / real_t(3.0) ;
	 const real_t CS_020 = CS_200 - CS_200__m__CS020 ;
	 const real_t CS_002 = CS_200 - CS_200__m__CS002 ;
	 
	 const real_t CS_120__p__CS_102 = (omega_trm3) * (C_120 + C_102) ;
	 const real_t CS_210__p__CS_012 = (omega_trm3) * (C_210 + C_012) ;
	 const real_t CS_201__p__CS_021 = (omega_trm3) * (C_201 + C_021) ;
	 const real_t CS_120__m__CS_102 = (omega_trm4) * (C_120 - C_102) ;
	 const real_t CS_210__m__CS_012 = (omega_trm4) * (C_210 - C_012) ;
	 const real_t CS_201__m__CS_021 = (omega_trm4) * (C_201 - C_021) ;
	 
	 const real_t CS_120 = real_t(0.5) * (CS_120__p__CS_102 + CS_120__m__CS_102) ;
	 const real_t CS_102 = real_t(0.5) * (CS_120__p__CS_102 - CS_120__m__CS_102) ;
	 const real_t CS_012 = real_t(0.5) * (CS_210__p__CS_012 - CS_210__m__CS_012) ;
	 const real_t CS_210 = real_t(0.5) * (CS_210__p__CS_012 + CS_210__m__CS_012) ;
	 const real_t CS_201 = real_t(0.5) * (CS_201__p__CS_021 + CS_201__m__CS_021) ;
	 const real_t CS_021 = real_t(0.5) * (CS_201__p__CS_021 - CS_201__m__CS_021) ;
	 const real_t CS_111 = (omega_trm5) * C_111 ;
	 
	 const real_t CS_220__m__2CS_202__p__CS_022 = (omega_trm6) * (C_220 - real_t(2.0) * C_202 + C_022) ;
	 const real_t CS_220__p__CS_202__m__2CS_022 = (omega_trm6) * (C_220 + C_202 - real_t(2.0) * C_022) ; 
	 const real_t CS_220__p__CS_202__p__CS_022  = (omega_trm7) * (C_220 + C_202 + C_022) ;
	 
	 const real_t CS_220 = (CS_220__m__2CS_202__p__CS_022 + CS_220__p__CS_202__m__2CS_022 + CS_220__p__CS_202__p__CS_022) / real_t(3.0) ;
	 const real_t CS_202 = (CS_220__p__CS_202__p__CS_022 - CS_220__m__2CS_202__p__CS_022) / real_t(3.0) ;
	 const real_t CS_022 = (CS_220__p__CS_202__p__CS_022 - CS_220__p__CS_202__m__2CS_022) / real_t(3.0) ;
	 
	 const real_t CS_211 = (omega_trm8 ) * C_211 ;
	 const real_t CS_121 = (omega_trm8 ) * C_121 ;
	 const real_t CS_112 = (omega_trm8) * C_112 ;
	 const real_t CS_221 = (omega_trm9 ) * C_221 ;
	 const real_t CS_212 = (omega_trm9 ) * C_212 ;
	 const real_t CS_122 = (omega_trm9 ) * C_122 ;
	 const real_t CS_222 = (omega_trm10) * C_222 ;
	 
	 // transforming after collision 
 
	 const real_t KC_000 = CS_000;
	 const real_t KC_100 = CS_100;
	 const real_t KC_200 = CS_200;
	 const real_t KC_001 = CS_001;
	 const real_t KC_101 = CS_101;
	 const real_t KC_201 = CS_201;
	 const real_t KC_002 = CS_002;
	 const real_t KC_102 = CS_102;
	 
	 const real_t sqr_KC_101 = KC_101 * KC_101 ;
	 const real_t KC_202 = CS_202 + (KC_002 * KC_200 + real_t(2.0) * sqr_KC_101) * rho_inv ;
	 const real_t KC_010 = CS_010 ;
	 const real_t KC_110 = CS_110 ;
	 const real_t KC_210 = CS_210 ;
	 const real_t KC_011 = CS_011 ;
	 const real_t KC_111 = CS_111 ;
	 const real_t KC_211 = CS_211 + (KC_200 * KC_011 + real_t(2.0) * KC_110 * KC_101) * rho_inv ;
	 const real_t KC_012 = CS_012;
	 const real_t KC_112 = CS_112 + (KC_002 * KC_110 + real_t(2.0) * KC_101 * KC_011) * rho_inv ;
	 const real_t KC_212 = CS_212 + (KC_002 * KC_210 + KC_200 * KC_210 + real_t(4.0) * KC_101 * KC_111 + real_t(2.0) * (KC_011 * KC_201 + KC_110 * KC_102)) * rho_inv;
	 const real_t KC_020 = CS_020;
	 const real_t KC_120 = CS_120;
	 
	 const real_t sqr_KC_110 = KC_110 * KC_110 ;
	 const real_t KC_220 = CS_220 + (KC_200 * KC_020 + real_t(2.0) * sqr_KC_110) * rho_inv;
	 const real_t KC_021 = CS_021;
	 const real_t KC_121 = CS_121 + (KC_020 * KC_101 + real_t(2.0) * KC_011 * KC_110) * rho_inv;
	 const real_t KC_221 = CS_221 + (KC_020 * KC_201 + KC_200 * KC_201 + real_t(4.0) * KC_110 * KC_111 + real_t(2.0) * (KC_011 * KC_210 + KC_101 * KC_120)) * rho_inv;
	 
	 const real_t sqr_KC_011 = KC_011 * KC_011 ;
	 const real_t KC_022 = CS_022 + (KC_020 * KC_002 + real_t(2.0) * sqr_KC_011) * rho_inv;
	 const real_t KC_122 = CS_122 + (KC_002 * KC_120 + KC_020 * KC_120 + real_t(4.0) * KC_011 * KC_111 + real_t(2.0) * (KC_101 * KC_021 + KC_110 * KC_012)) * rho_inv;
	 
	 const real_t sqr_KC_111 = KC_111 * KC_111 ;
	 const real_t KC_222 = CS_222 + (real_t(4.0) * sqr_KC_111 + KC_200 * KC_022 + KC_020 * KC_202 + KC_002 * KC_220 + real_t(4.0) * (KC_011 * KC_211 + KC_101 * KC_121 + KC_110 * KC_112) + real_t(2.0) * (KC_120 * KC_102 + KC_210 * KC_012 + KC_201 * KC_021)) * rho_inv
	 - (real_t(16.0) * KC_110 * KC_101 * KC_011 + real_t(4.0) * (KC_020 * sqr_KC_101 + KC_200 * sqr_KC_011 + KC_002 * sqr_KC_110) + real_t(2.0) * KC_200 * KC_020 * KC_002) * sqr_rho_inv;
                     
	 // trnasform back to central moment space 
	 // transform from central moment space to distribution funtion 
	 
	 // const defined for velocity in X direction 
	 const real_t oneminus_sqr_velx = real_t(1.0) - velXX ;
	 const real_t sqr_velx_plus_velX = velXX + velX ;
	 const real_t sqr_velx_minus_velX = velXX - velX ;
	 const real_t velx_term_plus = real_t(2.0) * velX + real_t(1.0) ;
	 const real_t velx_term_minus = real_t(2.0) * velX - real_t(1.0) ;
	 
	 
	 const real_t KC_1_00 = KC_000 * oneminus_sqr_velx	- KC_100 * real_t(2.0) * velX		- KC_200;
	 const real_t KC_0_00 = (KC_000 * sqr_velx_minus_velX	+ KC_100 * velx_term_minus		+ KC_200) * real_t(0.5);
	 const real_t KC_2_00 = (KC_000 * sqr_velx_plus_velX	+ KC_100 * velx_term_plus		+ KC_200) * real_t(0.5);
	 const real_t KC_1_01 = KC_001 * oneminus_sqr_velx	- KC_101 * real_t(2.0) * velX		- KC_201;
	 const real_t KC_0_01 = (KC_001 * sqr_velx_minus_velX	+ KC_101 * velx_term_minus		+ KC_201) * real_t(0.5);
	 const real_t KC_2_01 = (KC_001 * sqr_velx_plus_velX	+ KC_101 * velx_term_plus		+ KC_201) * real_t(0.5);
	 const real_t KC_1_02 = KC_002 * oneminus_sqr_velx	- KC_102 * real_t(2.0) * velX		- KC_202;
	 const real_t KC_0_02 = (KC_002 * sqr_velx_minus_velX	+ KC_102 * velx_term_minus		+ KC_202) * real_t(0.5);
	 const real_t KC_2_02 = (KC_002 * sqr_velx_plus_velX	+ KC_102 * velx_term_plus		+ KC_202) * real_t(0.5);
	 const real_t KC_1_10 = KC_010 * oneminus_sqr_velx	- KC_110 * real_t(2.0) * velX		- KC_210;
	 const real_t KC_0_10 = (KC_010 * sqr_velx_minus_velX	+ KC_110 * velx_term_minus		+ KC_210) * real_t(0.5);
	 const real_t KC_2_10 = (KC_010 * sqr_velx_plus_velX	+ KC_110 * velx_term_plus		+ KC_210) * real_t(0.5);
	 const real_t KC_1_11 = KC_011 * oneminus_sqr_velx	- KC_111 * real_t(2.0) * velX		- KC_211;
	 const real_t KC_0_11 = (KC_011 * sqr_velx_minus_velX	+ KC_111 * velx_term_minus		+ KC_211) * real_t(0.5);
	 const real_t KC_2_11 = (KC_011 * sqr_velx_plus_velX	+ KC_111 * velx_term_plus		+ KC_211) * real_t(0.5);
	 const real_t KC_1_12 = KC_012 * oneminus_sqr_velx	- KC_112 * real_t(2.0) * velX		- KC_212;
	 const real_t KC_0_12 = (KC_012 * sqr_velx_minus_velX	+ KC_112 * velx_term_minus		+ KC_212) * real_t(0.5);
	 const real_t KC_2_12 = (KC_012 * sqr_velx_plus_velX	+ KC_112 * velx_term_plus		+ KC_212) * real_t(0.5);
	 const real_t KC_1_20 = KC_020 * oneminus_sqr_velx	- KC_120 * real_t(2.0) * velX		- KC_220;
	 const real_t KC_0_20 = (KC_020 * sqr_velx_minus_velX	+ KC_120 * velx_term_minus		+ KC_220) * real_t(0.5);
	 const real_t KC_2_20 = (KC_020 * sqr_velx_plus_velX	+ KC_120 * velx_term_plus		+ KC_220) * real_t(0.5);
	 const real_t KC_1_21 = KC_021 * oneminus_sqr_velx	- KC_121 * real_t(2.0) * velX		- KC_221;
	 const real_t KC_0_21 = (KC_021 * sqr_velx_minus_velX	+ KC_121 * velx_term_minus		+ KC_221) * real_t(0.5);
	 const real_t KC_2_21 = (KC_021 * sqr_velx_plus_velX	+ KC_121 * velx_term_plus		+ KC_221) * real_t(0.5);
	 const real_t KC_1_22 = KC_022 * oneminus_sqr_velx	- KC_122 * real_t(2.0) * velX		- KC_222;
	 const real_t KC_0_22 = (KC_022 * sqr_velx_minus_velX	+ KC_122 * velx_term_minus		+ KC_222) * real_t(0.5);
	 const real_t KC_2_22 = (KC_022 * sqr_velx_plus_velX	+ KC_122 * velx_term_plus		+ KC_222) * real_t(0.5);
	
	
	 // collision is taking place from here and I need to change it from here 	
         // transform from velocity space to moment space and then to cumulant space , perform collsiopn and then again back transform to velocity space // 
         
         // const defined for velocity in Y direction 
         const real_t oneminus_sqr_vely = real_t(1.0) - velYY ;
	 const real_t sqr_vely_plus_velY = velYY + velY ;
	 const real_t sqr_vely_minus_velY = velYY - velY ;
	 const real_t vely_term_plus = real_t(2.0) * velY + real_t(1.0) ;
	 const real_t vely_term_minus = real_t(2.0) * velY - real_t(1.0) ;
	 
	 
	 const real_t KC_01_0 = KC_0_00 * oneminus_sqr_vely	- KC_0_10 * real_t(2.0) * velY		- KC_0_20;
	 const real_t KC_00_0 = (KC_0_00 * sqr_vely_minus_velY	+ KC_0_10 * vely_term_minus		+ KC_0_20) * real_t(0.5);
	 const real_t KC_02_0 = (KC_0_00 * sqr_vely_plus_velY	+ KC_0_10 * vely_term_plus		+ KC_0_20) * real_t(0.5);
	 const real_t KC_01_1 = KC_0_01 * oneminus_sqr_vely	- KC_0_11 * real_t(2.0) * velY		- KC_0_21;
	 const real_t KC_00_1 = (KC_0_01 * sqr_vely_minus_velY	+ KC_0_11 * vely_term_minus		+ KC_0_21) * real_t(0.5);
	 const real_t KC_02_1 = (KC_0_01 * sqr_vely_plus_velY	+ KC_0_11 * vely_term_plus		+ KC_0_21) * real_t(0.5);
	 const real_t KC_01_2 = KC_0_02 * oneminus_sqr_vely	- KC_0_12 * real_t(2.0) * velY		- KC_0_22;
	 const real_t KC_00_2 = (KC_0_02 * sqr_vely_minus_velY	+ KC_0_12 * vely_term_minus		+ KC_0_22) * real_t(0.5);
	 const real_t KC_02_2 = (KC_0_02 * sqr_vely_plus_velY	+ KC_0_12 * vely_term_plus		+ KC_0_22) * real_t(0.5);
	 const real_t KC_11_0 = KC_1_00 * oneminus_sqr_vely	- KC_1_10 * real_t(2.0) * velY		- KC_1_20;
	 const real_t KC_10_0 = (KC_1_00 * sqr_vely_minus_velY	+ KC_1_10 * vely_term_minus		+ KC_1_20) * real_t(0.5);
	 const real_t KC_12_0 = (KC_1_00 * sqr_vely_plus_velY	+ KC_1_10 * vely_term_plus		+ KC_1_20) * real_t(0.5);
	 const real_t KC_11_1 = KC_1_01 * oneminus_sqr_vely	- KC_1_11 * real_t(2.0) * velY		- KC_1_21;
	 const real_t KC_10_1 = (KC_1_01 * sqr_vely_minus_velY	+ KC_1_11 * vely_term_minus		+ KC_1_21) * real_t(0.5);
	 const real_t KC_12_1 = (KC_1_01 * sqr_vely_plus_velY	+ KC_1_11 * vely_term_plus		+ KC_1_21) * real_t(0.5);
	 const real_t KC_11_2 = KC_1_02 * oneminus_sqr_vely	- KC_1_12 * real_t(2.0) * velY		- KC_1_22;
	 const real_t KC_10_2 = (KC_1_02 * sqr_vely_minus_velY	+ KC_1_12 * vely_term_minus		+ KC_1_22) * real_t(0.5);
	 const real_t KC_12_2 = (KC_1_02 * sqr_vely_plus_velY	+ KC_1_12 * vely_term_plus		+ KC_1_22) * real_t(0.5);
	 const real_t KC_21_0 = KC_2_00 * oneminus_sqr_vely	- KC_2_10 * real_t(2.0) * velY		- KC_2_20;
	 const real_t KC_20_0 = (KC_2_00 * sqr_vely_minus_velY	+ KC_2_10 * vely_term_minus		+ KC_2_20) * real_t(0.5);
	 const real_t KC_22_0 = (KC_2_00 * sqr_vely_plus_velY	+ KC_2_10 * vely_term_plus		+ KC_2_20) * real_t(0.5);
	 const real_t KC_21_1 = KC_2_01 * oneminus_sqr_vely	- KC_2_11 * real_t(2.0) * velY		- KC_2_21;
	 const real_t KC_20_1 = (KC_2_01 * sqr_vely_minus_velY	+ KC_2_11 * vely_term_minus		+ KC_2_21) * real_t(0.5);
	 const real_t KC_22_1 = (KC_2_01 * sqr_vely_plus_velY	+ KC_2_11 * vely_term_plus		+ KC_2_21) * real_t(0.5);
	 const real_t KC_21_2 = KC_2_02 * oneminus_sqr_vely	- KC_2_12 * real_t(2.0) * velY		- KC_2_22;
	 const real_t KC_20_2 = (KC_2_02 * sqr_vely_minus_velY	+ KC_2_12 * vely_term_minus		+ KC_2_22) * real_t(0.5);
	 const real_t KC_22_2 = (KC_2_02 * sqr_vely_plus_velY	+ KC_2_12 * vely_term_plus		+ KC_2_22) * real_t(0.5);


	 // const defined for velocity in Z direction 
	 const real_t oneminus_sqr_velz = real_t(1.0) - velZZ ;
	 const real_t sqr_velz_plus_velZ = velZZ + velZ ;
	 const real_t sqr_velz_minus_velZ = velZZ - velZ ;
	 const real_t velz_term_plus = real_t(2.0) * velZ + real_t(1.0) ;
	 const real_t velz_term_minus = real_t(2.0) * velZ - real_t(1.0) ;
	 
	 dst->get(x,y,z,Stencil_T::idx[SW])   = KC_00_0 * oneminus_sqr_velz	- KC_00_1 * real_t(2.0) * velZ		- KC_00_2;
	 dst->get(x,y,z,Stencil_T::idx[BSW])  = (KC_00_0 * sqr_velz_minus_velZ	+ KC_00_1 * velz_term_minus		+ KC_00_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TSW])  = (KC_00_0 * sqr_velz_plus_velZ	+ KC_00_1 * velz_term_plus		+ KC_00_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[W])    = KC_01_0 * oneminus_sqr_velz	- KC_01_1 * real_t(2.0) * velZ		- KC_01_2;
	 dst->get(x,y,z,Stencil_T::idx[BW])   = (KC_01_0 * sqr_velz_minus_velZ	+ KC_01_1 * velz_term_minus		+ KC_01_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TW])   = (KC_01_0 * sqr_velz_plus_velZ	+ KC_01_1 * velz_term_plus		+ KC_01_2) * real_t(0.5);	 
	 dst->get(x,y,z,Stencil_T::idx[NW])   = KC_02_0 * oneminus_sqr_velz	- KC_02_1 * real_t(2.0) * velZ		- KC_02_2;
	 dst->get(x,y,z,Stencil_T::idx[BNW])  = (KC_02_0 * sqr_velz_minus_velZ	+ KC_02_1 * velz_term_minus		+ KC_02_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TNW])  = (KC_02_0 * sqr_velz_plus_velZ	+ KC_02_1 * velz_term_plus		+ KC_02_2) * real_t(0.5);	 
	 dst->get(x,y,z,Stencil_T::idx[S])    = KC_10_0 * oneminus_sqr_velz	- KC_10_1 * real_t(2.0) * velZ		- KC_10_2;
	 dst->get(x,y,z,Stencil_T::idx[BS])   = (KC_10_0 * sqr_velz_minus_velZ	+ KC_10_1 * velz_term_minus		+ KC_10_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TS])   = (KC_10_0 * sqr_velz_plus_velZ	+ KC_10_1 * velz_term_plus		+ KC_10_2) * real_t(0.5);	 	 
	 dst->get(x,y,z,Stencil_T::idx[C])    = KC_11_0 * oneminus_sqr_velz	- KC_11_1 * real_t(2.0) * velZ		- KC_11_2;
	 dst->get(x,y,z,Stencil_T::idx[B])    = (KC_11_0 * sqr_velz_minus_velZ	+ KC_11_1 * velz_term_minus		+ KC_11_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[T])    = (KC_11_0 * sqr_velz_plus_velZ	+ KC_11_1 * velz_term_plus		+ KC_11_2) * real_t(0.5);	 	 
	 dst->get(x,y,z,Stencil_T::idx[N])    = KC_12_0 * oneminus_sqr_velz	- KC_12_1 * real_t(2.0) * velZ		- KC_12_2;
	 dst->get(x,y,z,Stencil_T::idx[BN])   = (KC_12_0 * sqr_velz_minus_velZ	+ KC_12_1 * velz_term_minus		+ KC_12_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TN])   = (KC_12_0 * sqr_velz_plus_velZ	+ KC_12_1 * velz_term_plus		+ KC_12_2) * real_t(0.5);	 
	 dst->get(x,y,z,Stencil_T::idx[SE])   = KC_20_0 * oneminus_sqr_velz	- KC_20_1 * real_t(2.0) * velZ		- KC_20_2;
	 dst->get(x,y,z,Stencil_T::idx[BSE])  = (KC_20_0 * sqr_velz_minus_velZ	+ KC_20_1 * velz_term_minus		+ KC_20_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TSE])  = (KC_20_0 * sqr_velz_plus_velZ	+ KC_20_1 * velz_term_plus		+ KC_20_2) * real_t(0.5);	 
	 dst->get(x,y,z,Stencil_T::idx[E])    = KC_21_0 * oneminus_sqr_velz	- KC_21_1 * real_t(2.0) * velZ		- KC_21_2;
	 dst->get(x,y,z,Stencil_T::idx[BE])   = (KC_21_0 * sqr_velz_minus_velZ	+ KC_21_1 * velz_term_minus		+ KC_21_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TE])   = (KC_21_0 * sqr_velz_plus_velZ	+ KC_21_1 * velz_term_plus		+ KC_21_2) * real_t(0.5);	 
	 dst->get(x,y,z,Stencil_T::idx[NE])   = KC_22_0 * oneminus_sqr_velz	- KC_22_1 * real_t(2.0) * velZ		- KC_22_2;
	 dst->get(x,y,z,Stencil_T::idx[BNE])  = (KC_22_0 * sqr_velz_minus_velZ	+ KC_22_1 * velz_term_minus		+ KC_22_2) * real_t(0.5);
	 dst->get(x,y,z,Stencil_T::idx[TNE])  = (KC_22_0 * sqr_velz_plus_velZ	+ KC_22_1 * velz_term_plus		+ KC_22_2) * real_t(0.5);	 
	 
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
   
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_CUMULANT_1 )
{
      //WALBERLA_LOG_DEVEL("Cumulant Collide without streaming");
     
      // retrieving all omegas
      const real_t omega1 = src->latticeModel().collisionModel().omega1(); 
      const real_t omega2 = src->latticeModel().collisionModel().omega2();
      const real_t omega3 = src->latticeModel().collisionModel().omega3();
      const real_t omega4 = src->latticeModel().collisionModel().omega4();
      const real_t omega5 = src->latticeModel().collisionModel().omega5();
      const real_t omega6 = src->latticeModel().collisionModel().omega6();
      const real_t omega7 = src->latticeModel().collisionModel().omega7();
      const real_t omega8 = src->latticeModel().collisionModel().omega8();
      const real_t omega9 = src->latticeModel().collisionModel().omega9();
      const real_t omega10 = src->latticeModel().collisionModel().omega10();
   
      // This is the omega term used in the collision process 
      const real_t omega_trm1( real_t(1.0) - omega1 );
      const real_t omega_trm2( real_t(1.0) - omega2 );
      const real_t omega_trm3( real_t(1.0) - omega3 );
      const real_t omega_trm4( real_t(1.0) - omega4 );
      const real_t omega_trm5( real_t(1.0) - omega5 );
      const real_t omega_trm6( real_t(1.0) - omega6 );
      const real_t omega_trm7( real_t(1.0) - omega7 );
      const real_t omega_trm8( real_t(1.0) - omega8 );
      const real_t omega_trm9( real_t(1.0) - omega9 );
      const real_t omega_trm10( real_t(1.0) - omega10 );
  
  
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( this->filter(x,y,z) )
      {
         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP()
	 
	 const real_t updated_rho = rho;
	 
	 this->densityVelocityOut( x, y, z, lm, Vector3<real_t>( velX, velY, velZ ), updated_rho ); 
 
	 const real_t velXX = velX * velX;
         const real_t velYY = velY * velY;
         const real_t velZZ = velZ * velZ;
	 
	 // defining velocity terms for central moment space 
	 // velocity in z direction 
	 
	 const real_t velz_term = real_t(-1.0) * velZ ;
	 const real_t sqr_velz_term = velz_term * velz_term ;
	 const real_t oneminus_velz = real_t(1.0) - velZ ;
	 const real_t sqr_oneminus_velz = oneminus_velz * oneminus_velz ;
	 const real_t negoneminus_velz = real_t(-1.0) - velZ ;
	 const real_t sqr_negoneminus_velz = negoneminus_velz * negoneminus_velz ; 
	 
	 // Declaring the constants and initialising them with distribution function to transform  for the moment space 
	 const real_t k_00_0 = vBSW                        + vSW                    + vTSW ;
	 const real_t k_00_1 = vBSW * negoneminus_velz     + vSW * velz_term        + vTSW * oneminus_velz ;
	 const real_t k_00_2 = vBSW * sqr_negoneminus_velz + vSW * sqr_velz_term    + vTSW * sqr_oneminus_velz ;
	 const real_t k_01_0 = vBW                         + vW                     + vTW  ; 
	 const real_t k_01_1 = vBW * negoneminus_velz      + vW  * velz_term        + vTW * oneminus_velz ;
	 const real_t k_01_2 = vBW * sqr_negoneminus_velz  + vW  * sqr_velz_term    + vTW * sqr_oneminus_velz ;
	 const real_t k_02_0 = vBNW                        + vNW                    + vTNW ;
	 const real_t k_02_1 = vBNW * negoneminus_velz     + vNW * velz_term        + vTNW * oneminus_velz ;
	 const real_t k_02_2 = vBNW * sqr_negoneminus_velz + vNW * sqr_velz_term    + vTNW * sqr_oneminus_velz;
	 const real_t k_10_0 = vBS                         + vS                     + vTS;
	 const real_t k_10_1 = vBS * negoneminus_velz      + vS * velz_term         + vTS * oneminus_velz ;
	 const real_t k_10_2 = vBS * sqr_negoneminus_velz  + vS * sqr_velz_term     + vTS * sqr_oneminus_velz ;
	 const real_t k_11_0 = vB                          + vC                     + vT ;
	 const real_t k_11_1 = vB * negoneminus_velz       + vC * velz_term         + vT * oneminus_velz ;
	 const real_t k_11_2 = vB * sqr_negoneminus_velz   + vC * sqr_velz_term     + vT * sqr_oneminus_velz ;
	 const real_t k_12_0 = vBN                         + vN                     + vTN ;
	 const real_t k_12_1 = vBN * negoneminus_velz      + vN * velz_term         + vTN * oneminus_velz ;
	 const real_t k_12_2 = vBN * sqr_negoneminus_velz  + vN * sqr_velz_term     + vTN * sqr_oneminus_velz ;
	 const real_t k_20_0 = vBSE                        + vSE                    + vTSE ;
	 const real_t k_20_1 = vBSE * negoneminus_velz     + vSE * velz_term        + vTSE * oneminus_velz ;
	 const real_t k_20_2 = vBSE * sqr_negoneminus_velz + vSE * sqr_velz_term    + vTSE * sqr_oneminus_velz ;
	 const real_t k_21_0 = vBE                         + vE                     + vTE ;
	 const real_t k_21_1 = vBE * negoneminus_velz      + vE * velz_term         + vTE * oneminus_velz ;
	 const real_t k_21_2 = vBE * sqr_negoneminus_velz  + vE * sqr_velz_term     + vTE * sqr_oneminus_velz ;
	 const real_t k_22_0 = vBNE                        + vNE                    + vTNE ;
	 const real_t k_22_1 = vBNE * negoneminus_velz     + vNE * velz_term        + vTNE * oneminus_velz ;
	 const real_t k_22_2 = vBNE * sqr_negoneminus_velz + vNE * sqr_velz_term    + vTNE * sqr_oneminus_velz ;
	 
	  // velocity in y direction 
	 const real_t vely_term = real_t(-1.0) * velY ;
	 const real_t sqr_vely_term = vely_term * vely_term ;
	 const real_t oneminus_vely = real_t(1.0) - velY ;
	 const real_t sqr_oneminus_vely = oneminus_vely * oneminus_vely ;
	 const real_t negoneminus_vely = real_t(-1.0) - velY ;
	 const real_t sqr_negoneminus_vely = negoneminus_vely * negoneminus_vely ;
	 
	 const real_t k_0_00 = k_00_0				+ k_01_0			+ k_02_0;
	 const real_t k_0_10 = k_00_0 * negoneminus_vely	+ k_01_0 * vely_term		+ k_02_0 * oneminus_vely;
	 const real_t k_0_20 = k_00_0 * sqr_negoneminus_vely	+ k_01_0 * sqr_vely_term	+ k_02_0 * sqr_oneminus_vely; 
	 const real_t k_0_01 = k_00_1				+ k_01_1			+ k_02_1;
	 const real_t k_0_11 = k_00_1 * negoneminus_vely	+ k_01_1 * vely_term		+ k_02_1 *  oneminus_vely;
	 const real_t k_0_21 = k_00_1 * sqr_negoneminus_vely	+ k_01_1 * sqr_vely_term	+ k_02_1 * sqr_oneminus_vely;
	 const real_t k_0_02 = k_00_2				+ k_01_2			+ k_02_2;
	 const real_t k_0_12 = k_00_2 * negoneminus_vely	+ k_01_2 * vely_term		+ k_02_2 *  oneminus_vely;
	 const real_t k_0_22 = k_00_2 * sqr_negoneminus_vely	+ k_01_2 * sqr_vely_term	+ k_02_2 * sqr_oneminus_vely;
	 const real_t k_1_00 = k_10_0				+ k_11_0			+ k_12_0;
	 const real_t k_1_10 = k_10_0 * negoneminus_vely	+ k_11_0 * vely_term		+ k_12_0 *  oneminus_vely;
	 const real_t k_1_20 = k_10_0 * sqr_negoneminus_vely	+ k_11_0 * sqr_vely_term	+ k_12_0 * sqr_oneminus_vely;
	 const real_t k_1_01 = k_10_1				+ k_11_1			+ k_12_1;
	 const real_t k_1_11 = k_10_1 * negoneminus_vely	+ k_11_1 * vely_term		+ k_12_1 *  oneminus_vely;
	 const real_t k_1_21 = k_10_1 * sqr_negoneminus_vely	+ k_11_1 * sqr_vely_term	+ k_12_1 * sqr_oneminus_vely;
	 const real_t k_1_02 = k_10_2 				+ k_11_2			+ k_12_2 ;
	 const real_t k_1_12 = k_10_2 * negoneminus_vely	+ k_11_2 * vely_term		+ k_12_2 *  oneminus_vely;
	 const real_t k_1_22 = k_10_2 * sqr_negoneminus_vely	+ k_11_2 * sqr_vely_term	+ k_12_2 * sqr_oneminus_vely;
	 const real_t k_2_00 = k_20_0				+ k_21_0			+ k_22_0;
	 const real_t k_2_10 = k_20_0 * negoneminus_vely	+ k_21_0 * vely_term		+ k_22_0 *  oneminus_vely;
	 const real_t k_2_20 = k_20_0 * sqr_negoneminus_vely	+ k_21_0 * sqr_vely_term	+ k_22_0 * sqr_oneminus_vely;
	 const real_t k_2_01 = k_20_1				+ k_21_1			+ k_22_1;
	 const real_t k_2_11 = k_20_1 * negoneminus_vely	+ k_21_1 * vely_term		+ k_22_1 *  oneminus_vely;
	 const real_t k_2_21 = k_20_1 * sqr_negoneminus_vely	+ k_21_1 * sqr_vely_term	+ k_22_1 * sqr_oneminus_vely;
	 const real_t k_2_02 = k_20_2				+ k_21_2			+ k_22_2;
	 const real_t k_2_12 = k_20_2 * negoneminus_vely	+ k_21_2 * vely_term		+ k_22_2 *  oneminus_vely;
	 const real_t k_2_22 = k_20_2 * sqr_negoneminus_vely	+ k_21_2 * sqr_vely_term	+ k_22_2 * sqr_oneminus_vely;
	 
	  // Velocity in x direction 
	 const real_t velx_term = real_t(-1.0) * velX ;
	 const real_t sqr_velx_term = velx_term * velx_term ;
	 const real_t oneminus_velx = real_t(1.0) - velX ;
	 const real_t sqr_oneminus_velx = oneminus_velx * oneminus_velx ;
	 const real_t negoneminus_velx = real_t(-1.0) - velX ;
	 const real_t sqr_negoneminus_velx = negoneminus_velx * negoneminus_velx ;
	 
	 const real_t k_000 = k_0_00				+ k_1_00			+ k_2_00 ;
	 const real_t k_100 = k_0_00 * negoneminus_velx		+ k_1_00 * velx_term		+ k_2_00 * oneminus_velx ;
	 const real_t k_200 = k_0_00 * sqr_negoneminus_velx	+ k_1_00 * sqr_velx_term	+ k_2_00 * sqr_oneminus_velx ;
	 const real_t k_001 = k_0_01				+ k_1_01			+ k_2_01 ;
	 const real_t k_101 = k_0_01 * negoneminus_velx		+ k_1_01 * velx_term		+ k_2_01 * oneminus_velx;
	 const real_t k_201 = k_0_01 * sqr_negoneminus_velx	+ k_1_01 * sqr_velx_term	+ k_2_01 * sqr_oneminus_velx;
	 const real_t k_002 = k_0_02				+ k_1_02			+ k_2_02;
	 const real_t k_102 = k_0_02 * negoneminus_velx		+ k_1_02 * velx_term		+ k_2_02 * oneminus_velx;
	 const real_t k_202 = k_0_02 * sqr_negoneminus_velx	+ k_1_02 * sqr_velx_term	+ k_2_02 * sqr_oneminus_velx;
	 const real_t k_010 = k_0_10				+ k_1_10			+ k_2_10;
	 const real_t k_110 = k_0_10 * negoneminus_velx		+ k_1_10 * velx_term		+ k_2_10 * oneminus_velx;
	 const real_t k_210 = k_0_10 * sqr_negoneminus_velx	+ k_1_10 * sqr_velx_term	+ k_2_10 * sqr_oneminus_velx;
	 const real_t k_011 = k_0_11				+ k_1_11			+ k_2_11;
	 const real_t k_111 = k_0_11 * negoneminus_velx		+ k_1_11 * velx_term		+ k_2_11 * oneminus_velx;
	 const real_t k_211 = k_0_11 * sqr_negoneminus_velx	+ k_1_11 * sqr_velx_term	+ k_2_11 * sqr_oneminus_velx;
	 const real_t k_012 = k_0_12				+ k_1_12			+ k_2_12;
	 const real_t k_112 = k_0_12 * negoneminus_velx		+ k_1_12 * velx_term		+ k_2_12 * oneminus_velx;
	 const real_t k_212 = k_0_12 * sqr_negoneminus_velx	+ k_1_12 * sqr_velx_term	+ k_2_12 * sqr_oneminus_velx;
	 const real_t k_020 = k_0_20				+ k_1_20			+ k_2_20;
	 const real_t k_120 = k_0_20 * negoneminus_velx		+ k_1_20 * velx_term		+ k_2_20 * oneminus_velx;
	 const real_t k_220 = k_0_20 * sqr_negoneminus_velx	+ k_1_20 * sqr_velx_term	+ k_2_20 * sqr_oneminus_velx;
	 const real_t k_021 = k_0_21				+ k_1_21			+ k_2_21;
	 const real_t k_121 = k_0_21 * negoneminus_velx		+ k_1_21 * velx_term		+ k_2_21 * oneminus_velx;
	 const real_t k_221 = k_0_21 * sqr_negoneminus_velx	+ k_1_21 * sqr_velx_term	+ k_2_21 * sqr_oneminus_velx;
	 const real_t k_022 = k_0_22				+ k_1_22			+ k_2_22;
	 const real_t k_122 = k_0_22 * negoneminus_velx		+ k_1_22 * velx_term		+ k_2_22 * oneminus_velx;
	 const real_t k_222 = k_0_22 * sqr_negoneminus_velx	+ k_1_22 * sqr_velx_term	+ k_2_22 * sqr_oneminus_velx;
	 
	 // transformation to moment space done now , transform to cumulant space 
	 
	 // defining the constants for central moment space 
	 const real_t rho_inv = real_t(1.0) / updated_rho; 

	 // defining the sqaures 
	 const real_t sqr_k_110 = k_110 * k_110 ;
	 const real_t sqr_k_101 = k_101 * k_101 ;
	 const real_t sqr_k_011 = k_011 * k_011 ;
	 const real_t sqr_k_111 = k_111 * k_111 ;
	 
	 const real_t sqr_rho_inv = rho_inv * rho_inv ;
	 
	 // defining the cumulant space constants
	 const real_t C_000 = k_000 ;
	 const real_t C_001 = k_001 ;
	 const real_t C_002 = k_002 ;
	 const real_t C_010 = k_010 ;
	 const real_t C_011 = k_011 ;
	 const real_t C_012 = k_012 ;
	 const real_t C_020 = k_020 ;
	 const real_t C_021 = k_021 ;
	 const real_t C_022 = k_022 - (k_020 * k_002 + real_t(2.0) * sqr_k_011) * rho_inv ;
	 const real_t C_100 = k_100 ;
	 const real_t C_101 = k_101 ;
	 const real_t C_102 = k_102 ;
	 const real_t C_110 = k_110 ;
	 const real_t C_111 = k_111 ;
	 const real_t C_112 = k_112 - (k_002 * k_110 + real_t(2.0) * k_101 * k_011) * rho_inv ;
	 const real_t C_120 = k_120 ;
	 const real_t C_121 = k_121 - (k_020 * k_101 + real_t(2.0) * k_011 * k_110) * rho_inv ;
	 const real_t C_122 = k_122 - (k_002 * k_120 + k_020 * k_120 + real_t(4.0) * k_011 * k_111 + real_t(2.0) * (k_101 * k_021 + k_110 * k_012)) * rho_inv ;
	 const real_t C_200 = k_200 ;
	 const real_t C_201 = k_201 ;
	 const real_t C_202 = k_202 - (k_002 * k_200 + real_t(2.0) * sqr_k_101) * rho_inv ;
	 const real_t C_210 = k_210 ;
	 const real_t C_211 = k_211 - (k_200 * k_011 + real_t(2.0) * k_110 * k_101) * rho_inv ;
	 const real_t C_212 = k_212 - (k_002 * k_210 + k_200 * k_210 + real_t(4.0) * k_101 * k_111 + real_t(2.0) * (k_011 * k_201 + k_110 * k_102)) * rho_inv ;
	 const real_t C_220 = k_220 - (k_200 * k_020 + real_t(2.0) * sqr_k_110) * rho_inv ;
	 const real_t C_221 = k_221 - (k_020 * k_201 + k_200 * k_201 + real_t(4.0) * k_110 * k_111 + real_t(2.0) * (k_011 * k_210 + k_101 * k_120)) * rho_inv ;
	 const real_t C_222 = k_222 - (real_t(4.0)* sqr_k_111 + k_200 * k_022 + k_020 * k_202 + k_002 * k_220  + real_t(4.0) * (k_011 * k_211 + k_101 * k_121 + k_110 * k_112) + real_t(2.0) * (k_120 * k_102 + k_210 * k_012 + k_201 * k_021)) * rho_inv 
 				    + (real_t(16.0) * k_110 * k_101 * k_011 + real_t(4.0) * (k_020 * sqr_k_101 + k_200 * sqr_k_011 + k_002 * sqr_k_110) + real_t(2.0) * k_200 *k_020 * k_002) * sqr_rho_inv ;
	 
	 
	 // collision happens in cumulant space 
	 const real_t CS_000 = C_000 ;
	 const real_t CS_100 = C_100 ;
	 const real_t CS_010 = C_010 ;
	 const real_t CS_001 = C_001 ;
	 const real_t CS_110 = omega_trm1 * C_110 ;
	 const real_t CS_101 = omega_trm1 * C_101 ;
         const real_t CS_011 = omega_trm1 * C_011 ;
	 
	 const real_t Dxux = -real_t(0.5) * omega1 * rho_inv * (real_t(2.0) * C_200 - C_020 - C_002) - real_t(0.5) * omega2 * rho_inv * (C_200 + C_020 + C_002 - C_000) ;
         const real_t Dyuy = Dxux + real_t(1.5) * omega1 * rho_inv * (C_200 - C_020) ;
         const real_t Dzuz = Dxux + real_t(1.5) * omega1 * rho_inv * (C_200 - C_002) ;
	 
	 const real_t CS_200__m__CS020 = (omega_trm1) * (C_200 - C_020) - real_t(3.0) * updated_rho * (real_t(1.0) - real_t(0.5) * omega1) * (Dxux * velXX - Dyuy * velYY) ;
	 const real_t CS_200__m__CS002 = (omega_trm1) * (C_200 - C_002) - real_t(3.0) * updated_rho * (real_t(1.0) - real_t(0.5) * omega1) * (Dxux * velXX - Dzuz * velZZ) ;
	 const real_t CS_200__p__CS020__p__CS_002 = omega2 * C_000 + (omega_trm2) * (C_200 + C_020 + C_002) - real_t(3.0) * updated_rho *(real_t(1.0) - real_t(0.5) * omega2) * (Dxux * velXX + Dyuy * velYY + Dzuz * velZZ);
	 const real_t CS_200 = (CS_200__m__CS020 + CS_200__m__CS002 + CS_200__p__CS020__p__CS_002) / real_t(3.0) ;
         const real_t CS_020 = CS_200 - CS_200__m__CS020 ;
         const real_t CS_002 = CS_200 - CS_200__m__CS002 ;
	 const real_t CS_120__p__CS_102 = (omega_trm3) * (C_120 + C_102) ;
         const real_t CS_210__p__CS_012 = (omega_trm3) * (C_210 + C_012) ;
         const real_t CS_201__p__CS_021 = (omega_trm3) * (C_201 + C_021) ;
         const real_t CS_120__m__CS_102 = (omega_trm4) * (C_120 - C_102) ;
         const real_t CS_210__m__CS_012 = (omega_trm4) * (C_210 - C_012) ;
         const real_t CS_201__m__CS_021 = (omega_trm4) * (C_201 - C_021) ;
	 
	 const real_t CS_120 = real_t(0.5) * (CS_120__p__CS_102 + CS_120__m__CS_102) ;
	 const real_t CS_102 = real_t(0.5) * (CS_120__p__CS_102 - CS_120__m__CS_102) ;
         const real_t CS_012 = real_t(0.5) * (CS_210__p__CS_012 - CS_210__m__CS_012) ;
	 const real_t CS_210 = real_t(0.5) * (CS_210__p__CS_012 + CS_210__m__CS_012) ;
	 const real_t CS_201 = real_t(0.5) * (CS_201__p__CS_021 + CS_201__m__CS_021) ;
	 const real_t CS_021 = real_t(0.5) * (CS_201__p__CS_021 - CS_201__m__CS_021) ;
	 const real_t CS_111 = (omega_trm5) * C_111 ;
	 
	 const real_t CS_220__m__2CS_202__p__CS_022 = (omega_trm6) * (C_220 - real_t(2.0) * C_202 + C_022) ;
	 const real_t CS_220__p__CS_202__m__2CS_022 = (omega_trm6) * (C_220 + C_202 - real_t(2.0) * C_022) ; 
	 const real_t CS_220__p__CS_202__p__CS_022  = (omega_trm7) * (C_220 + C_202 + C_022) ;
	 
	 const real_t CS_220 = (CS_220__m__2CS_202__p__CS_022 + CS_220__p__CS_202__m__2CS_022 + CS_220__p__CS_202__p__CS_022) / real_t(3.0) ;
	 const real_t CS_202 = (CS_220__p__CS_202__p__CS_022 - CS_220__m__2CS_202__p__CS_022) / real_t(3.0) ;
	 const real_t CS_022 = (CS_220__p__CS_202__p__CS_022 - CS_220__p__CS_202__m__2CS_022) / real_t(3.0) ;
	 
	 const real_t CS_211 = (omega_trm8 ) * C_211 ;
	 const real_t CS_121 = (omega_trm8 ) * C_121 ;
	 const real_t CS_112 = (omega_trm8 ) * C_112 ;
	 const real_t CS_221 = (omega_trm9 ) * C_221 ;
	 const real_t CS_212 = (omega_trm9 ) * C_212 ;
	 const real_t CS_122 = (omega_trm9 ) * C_122 ;
	 const real_t CS_222 = (omega_trm10) * C_222 ;
	 
	 // transforming after collision 
 
	 const real_t KC_000 = CS_000;
	 const real_t KC_100 = CS_100;
	 const real_t KC_200 = CS_200;
	 const real_t KC_001 = CS_001;
	 const real_t KC_101 = CS_101;
	 const real_t KC_201 = CS_201;
	 const real_t KC_002 = CS_002;
	 const real_t KC_102 = CS_102;
	 
	 const real_t sqr_KC_101 = KC_101 * KC_101 ;
	 
	 const real_t KC_202 = CS_202 + (KC_002 * KC_200 + real_t(2.0) * sqr_KC_101) * rho_inv ;
	 const real_t KC_010 = CS_010 ;
	 const real_t KC_110 = CS_110 ;
	 const real_t KC_210 = CS_210 ;
	 const real_t KC_011 = CS_011 ;
	 const real_t KC_111 = CS_111 ;
	 const real_t KC_211 = CS_211 + (KC_200 * KC_011 + real_t(2.0) * KC_110 * KC_101) * rho_inv ;
	 const real_t KC_012 = CS_012;
	 const real_t KC_112 = CS_112 + (KC_002 * KC_110 + real_t(2.0) * KC_101 * KC_011) * rho_inv ;
	 const real_t KC_212 = CS_212 + (KC_002 * KC_210 + KC_200 * KC_210 + real_t(4.0) * KC_101 * KC_111 + real_t(2.0) * (KC_011 * KC_201 + KC_110 * KC_102)) * rho_inv ;
	 const real_t KC_020 = CS_020;
	 const real_t KC_120 = CS_120;
	 
	 const real_t sqr_KC_110 = KC_110 * KC_110 ;
	 
	 const real_t KC_220 = CS_220 + (KC_200 * KC_020 + real_t(2.0) * sqr_KC_110) * rho_inv;
	 const real_t KC_021 = CS_021;
	 const real_t KC_121 = CS_121 + (KC_020 * KC_101 + real_t(2.0) * KC_011 * KC_110) * rho_inv;
	 const real_t KC_221 = CS_221 + (KC_020 * KC_201 + KC_200 * KC_201 + real_t(4.0) * KC_110 * KC_111 + real_t(2.0) * (KC_011 * KC_210 + KC_101 * KC_120)) * rho_inv;
	 
	 const real_t sqr_KC_011 = KC_011 * KC_011 ;
	 
	 const real_t KC_022 = CS_022 + (KC_020 * KC_002 + real_t(2.0) * sqr_KC_011) * rho_inv;
	 const real_t KC_122 = CS_122 + (KC_002 * KC_120 + KC_020 * KC_120 + real_t(4.0) * KC_011 * KC_111 + real_t(2.0) * (KC_101 * KC_021 + KC_110 * KC_012)) * rho_inv;
	 
	 const real_t sqr_KC_111 = KC_111 * KC_111 ;
	 
	 const real_t KC_222 = CS_222 + (real_t(4.0) * sqr_KC_111 + KC_200 * KC_022 + KC_020 * KC_202 + KC_002 * KC_220 + real_t(4.0) * (KC_011 * KC_211 + KC_101 * KC_121 + KC_110 * KC_112) + real_t(2.0) * (KC_120 * KC_102 + KC_210 * KC_012 + KC_201 * KC_021)) * rho_inv
	 - (real_t(16.0) * KC_110 * KC_101 * KC_011 + real_t(4.0) * (KC_020 * sqr_KC_101 + KC_200 * sqr_KC_011 + KC_002 * sqr_KC_110) + real_t(2.0) * KC_200 * KC_020 * KC_002) * sqr_rho_inv;
                     
	 // trnasform back to central moment space 
	 // transform from central moment space to distribution funtion 
	 // const defined for velocity in X direction 
	 const real_t oneminus_sqr_velx = real_t(1.0) - velXX ;
	 const real_t sqr_velx_plus_velX = velXX + velX ;
	 const real_t sqr_velx_minus_velX = velXX - velX ;
	 const real_t velx_term_plus = real_t(2.0) * velX + real_t(1.0) ;
	 const real_t velx_term_minus = real_t(2.0) * velX - real_t(1.0) ;
	 
	 
	 const real_t KC_1_00 = KC_000 * oneminus_sqr_velx	- KC_100 * real_t(2.0) * velX		- KC_200;
	 const real_t KC_0_00 = (KC_000 * sqr_velx_minus_velX	+ KC_100 * velx_term_minus		+ KC_200) * real_t(0.5);
	 const real_t KC_2_00 = (KC_000 * sqr_velx_plus_velX	+ KC_100 * velx_term_plus		+ KC_200) * real_t(0.5);
	 const real_t KC_1_01 = KC_001 * oneminus_sqr_velx	- KC_101 * real_t(2.0) * velX		- KC_201;
	 const real_t KC_0_01 = (KC_001 * sqr_velx_minus_velX	+ KC_101 * velx_term_minus		+ KC_201) * real_t(0.5);
	 const real_t KC_2_01 = (KC_001 * sqr_velx_plus_velX	+ KC_101 * velx_term_plus		+ KC_201) * real_t(0.5);
	 const real_t KC_1_02 = KC_002 * oneminus_sqr_velx	- KC_102 * real_t(2.0) * velX		- KC_202;
	 const real_t KC_0_02 = (KC_002 * sqr_velx_minus_velX	+ KC_102 * velx_term_minus		+ KC_202) * real_t(0.5);
	 const real_t KC_2_02 = (KC_002 * sqr_velx_plus_velX	+ KC_102 * velx_term_plus		+ KC_202) * real_t(0.5);
	 const real_t KC_1_10 = KC_010 * oneminus_sqr_velx	- KC_110 * real_t(2.0) * velX		- KC_210;
	 const real_t KC_0_10 = (KC_010 * sqr_velx_minus_velX	+ KC_110 * velx_term_minus		+ KC_210) * real_t(0.5);
	 const real_t KC_2_10 = (KC_010 * sqr_velx_plus_velX	+ KC_110 * velx_term_plus		+ KC_210) * real_t(0.5);
	 const real_t KC_1_11 = KC_011 * oneminus_sqr_velx	- KC_111 * real_t(2.0) * velX		- KC_211;
	 const real_t KC_0_11 = (KC_011 * sqr_velx_minus_velX	+ KC_111 * velx_term_minus		+ KC_211) * real_t(0.5);
	 const real_t KC_2_11 = (KC_011 * sqr_velx_plus_velX	+ KC_111 * velx_term_plus		+ KC_211) * real_t(0.5);
	 const real_t KC_1_12 = KC_012 * oneminus_sqr_velx	- KC_112 * real_t(2.0) * velX		- KC_212;
	 const real_t KC_0_12 = (KC_012 * sqr_velx_minus_velX	+ KC_112 * velx_term_minus		+ KC_212) * real_t(0.5);
	 const real_t KC_2_12 = (KC_012 * sqr_velx_plus_velX	+ KC_112 * velx_term_plus		+ KC_212) * real_t(0.5);
	 const real_t KC_1_20 = KC_020 * oneminus_sqr_velx	- KC_120 * real_t(2.0) * velX		- KC_220;
	 const real_t KC_0_20 = (KC_020 * sqr_velx_minus_velX	+ KC_120 * velx_term_minus		+ KC_220) * real_t(0.5);
	 const real_t KC_2_20 = (KC_020 * sqr_velx_plus_velX	+ KC_120 * velx_term_plus		+ KC_220) * real_t(0.5);
	 const real_t KC_1_21 = KC_021 * oneminus_sqr_velx	- KC_121 * real_t(2.0) * velX		- KC_221;
	 const real_t KC_0_21 = (KC_021 * sqr_velx_minus_velX	+ KC_121 * velx_term_minus		+ KC_221) * real_t(0.5);
	 const real_t KC_2_21 = (KC_021 * sqr_velx_plus_velX	+ KC_121 * velx_term_plus		+ KC_221) * real_t(0.5);
	 const real_t KC_1_22 = KC_022 * oneminus_sqr_velx	- KC_122 * real_t(2.0) * velX		- KC_222;
	 const real_t KC_0_22 = (KC_022 * sqr_velx_minus_velX	+ KC_122 * velx_term_minus		+ KC_222) * real_t(0.5);
	 const real_t KC_2_22 = (KC_022 * sqr_velx_plus_velX	+ KC_122 * velx_term_plus		+ KC_222) * real_t(0.5);
	
	
	 // collision is taking place from here and I need to change it from here 	
         // transform from velocity space to moment space and then to cumulant space , perform collsiopn and then again back transform to velocity space // 
         
         // const defined for velocity in Y direction 
         const real_t oneminus_sqr_vely = real_t(1.0) - velYY ;
	 const real_t sqr_vely_plus_velY = velYY + velY ;
	 const real_t sqr_vely_minus_velY = velYY - velY ;
	 const real_t vely_term_plus = real_t(2.0) * velY + real_t(1.0) ;
	 const real_t vely_term_minus = real_t(2.0) * velY - real_t(1.0) ;
	 
	 
	 const real_t KC_01_0 = KC_0_00 * oneminus_sqr_vely	- KC_0_10 * real_t(2.0) * velY		- KC_0_20;
	 const real_t KC_00_0 = (KC_0_00 * sqr_vely_minus_velY	+ KC_0_10 * vely_term_minus		+ KC_0_20) * real_t(0.5);
	 const real_t KC_02_0 = (KC_0_00 * sqr_vely_plus_velY	+ KC_0_10 * vely_term_plus		+ KC_0_20) * real_t(0.5);
	 const real_t KC_01_1 = KC_0_01 * oneminus_sqr_vely	- KC_0_11 * real_t(2.0) * velY		- KC_0_21;
	 const real_t KC_00_1 = (KC_0_01 * sqr_vely_minus_velY	+ KC_0_11 * vely_term_minus		+ KC_0_21) * real_t(0.5);
	 const real_t KC_02_1 = (KC_0_01 * sqr_vely_plus_velY	+ KC_0_11 * vely_term_plus		+ KC_0_21) * real_t(0.5);
	 const real_t KC_01_2 = KC_0_02 * oneminus_sqr_vely	- KC_0_12 * real_t(2.0) * velY		- KC_0_22;
	 const real_t KC_00_2 = (KC_0_02 * sqr_vely_minus_velY	+ KC_0_12 * vely_term_minus		+ KC_0_22) * real_t(0.5);
	 const real_t KC_02_2 = (KC_0_02 * sqr_vely_plus_velY	+ KC_0_12 * vely_term_plus		+ KC_0_22) * real_t(0.5);
	 const real_t KC_11_0 = KC_1_00 * oneminus_sqr_vely	- KC_1_10 * real_t(2.0) * velY		- KC_1_20;
	 const real_t KC_10_0 = (KC_1_00 * sqr_vely_minus_velY	+ KC_1_10 * vely_term_minus		+ KC_1_20) * real_t(0.5);
	 const real_t KC_12_0 = (KC_1_00 * sqr_vely_plus_velY	+ KC_1_10 * vely_term_plus		+ KC_1_20) * real_t(0.5);
	 const real_t KC_11_1 = KC_1_01 * oneminus_sqr_vely	- KC_1_11 * real_t(2.0) * velY		- KC_1_21;
	 const real_t KC_10_1 = (KC_1_01 * sqr_vely_minus_velY	+ KC_1_11 * vely_term_minus		+ KC_1_21) * real_t(0.5);
	 const real_t KC_12_1 = (KC_1_01 * sqr_vely_plus_velY	+ KC_1_11 * vely_term_plus		+ KC_1_21) * real_t(0.5);
	 const real_t KC_11_2 = KC_1_02 * oneminus_sqr_vely	- KC_1_12 * real_t(2.0) * velY		- KC_1_22;
	 const real_t KC_10_2 = (KC_1_02 * sqr_vely_minus_velY	+ KC_1_12 * vely_term_minus		+ KC_1_22) * real_t(0.5);
	 const real_t KC_12_2 = (KC_1_02 * sqr_vely_plus_velY	+ KC_1_12 * vely_term_plus		+ KC_1_22) * real_t(0.5);
	 const real_t KC_21_0 = KC_2_00 * oneminus_sqr_vely	- KC_2_10 * real_t(2.0) * velY		- KC_2_20;
	 const real_t KC_20_0 = (KC_2_00 * sqr_vely_minus_velY	+ KC_2_10 * vely_term_minus		+ KC_2_20) * real_t(0.5);
	 const real_t KC_22_0 = (KC_2_00 * sqr_vely_plus_velY	+ KC_2_10 * vely_term_plus		+ KC_2_20) * real_t(0.5);
	 const real_t KC_21_1 = KC_2_01 * oneminus_sqr_vely	- KC_2_11 * real_t(2.0) * velY		- KC_2_21;
	 const real_t KC_20_1 = (KC_2_01 * sqr_vely_minus_velY	+ KC_2_11 * vely_term_minus		+ KC_2_21) * real_t(0.5);
	 const real_t KC_22_1 = (KC_2_01 * sqr_vely_plus_velY	+ KC_2_11 * vely_term_plus		+ KC_2_21) * real_t(0.5);
	 const real_t KC_21_2 = KC_2_02 * oneminus_sqr_vely	- KC_2_12 * real_t(2.0) * velY		- KC_2_22;
	 const real_t KC_20_2 = (KC_2_02 * sqr_vely_minus_velY	+ KC_2_12 * vely_term_minus		+ KC_2_22) * real_t(0.5);
	 const real_t KC_22_2 = (KC_2_02 * sqr_vely_plus_velY	+ KC_2_12 * vely_term_plus		+ KC_2_22) * real_t(0.5);


	 // const defined for velocity in Z direction 
	 const real_t oneminus_sqr_velz = real_t(1.0) - velZZ ;
	 const real_t sqr_velz_plus_velZ = velZZ + velZ ;
	 const real_t sqr_velz_minus_velZ = velZZ - velZ ;
	 const real_t velz_term_plus = real_t(2.0) * velZ + real_t(1.0) ;
	 const real_t velz_term_minus = real_t(2.0) * velZ - real_t(1.0) ;
	 
	 
	 // updating the distribution function 
	 src->get(x,y,z,Stencil_T::idx[SW])   = KC_00_0 * oneminus_sqr_velz	- KC_00_1 * real_t(2.0) * velZ		- KC_00_2;
	 src->get(x,y,z,Stencil_T::idx[BSW])  = (KC_00_0 * sqr_velz_minus_velZ	+ KC_00_1 * velz_term_minus		+ KC_00_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TSW])  = (KC_00_0 * sqr_velz_plus_velZ	+ KC_00_1 * velz_term_plus		+ KC_00_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[W])    = KC_01_0 * oneminus_sqr_velz	- KC_01_1 * real_t(2.0) * velZ		- KC_01_2;
	 src->get(x,y,z,Stencil_T::idx[BW])   = (KC_01_0 * sqr_velz_minus_velZ	+ KC_01_1 * velz_term_minus		+ KC_01_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TW])   = (KC_01_0 * sqr_velz_plus_velZ	+ KC_01_1 * velz_term_plus		+ KC_01_2) * real_t(0.5);	 
	 src->get(x,y,z,Stencil_T::idx[NW])   = KC_02_0 * oneminus_sqr_velz	- KC_02_1 * real_t(2.0) * velZ		- KC_02_2;
	 src->get(x,y,z,Stencil_T::idx[BNW])  = (KC_02_0 * sqr_velz_minus_velZ	+ KC_02_1 * velz_term_minus		+ KC_02_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TNW])  = (KC_02_0 * sqr_velz_plus_velZ	+ KC_02_1 * velz_term_plus		+ KC_02_2) * real_t(0.5);	 
	 src->get(x,y,z,Stencil_T::idx[S])    = KC_10_0 * oneminus_sqr_velz	- KC_10_1 * real_t(2.0) * velZ		- KC_10_2;
	 src->get(x,y,z,Stencil_T::idx[BS])   = (KC_10_0 * sqr_velz_minus_velZ	+ KC_10_1 * velz_term_minus		+ KC_10_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TS])   = (KC_10_0 * sqr_velz_plus_velZ	+ KC_10_1 * velz_term_plus		+ KC_10_2) * real_t(0.5);	 	 
	 src->get(x,y,z,Stencil_T::idx[C])    = KC_11_0 * oneminus_sqr_velz	- KC_11_1 * real_t(2.0) * velZ		- KC_11_2;
	 src->get(x,y,z,Stencil_T::idx[B])    = (KC_11_0 * sqr_velz_minus_velZ	+ KC_11_1 * velz_term_minus		+ KC_11_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[T])    = (KC_11_0 * sqr_velz_plus_velZ	+ KC_11_1 * velz_term_plus		+ KC_11_2) * real_t(0.5);	 	 
	 src->get(x,y,z,Stencil_T::idx[N])    = KC_12_0 * oneminus_sqr_velz	- KC_12_1 * real_t(2.0) * velZ		- KC_12_2;
	 src->get(x,y,z,Stencil_T::idx[BN])   = (KC_12_0 * sqr_velz_minus_velZ	+ KC_12_1 * velz_term_minus		+ KC_12_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TN])   = (KC_12_0 * sqr_velz_plus_velZ	+ KC_12_1 * velz_term_plus		+ KC_12_2) * real_t(0.5);	 
	 src->get(x,y,z,Stencil_T::idx[SE])   = KC_20_0 * oneminus_sqr_velz	- KC_20_1 * real_t(2.0) * velZ		- KC_20_2;
	 src->get(x,y,z,Stencil_T::idx[BSE])  = (KC_20_0 * sqr_velz_minus_velZ	+ KC_20_1 * velz_term_minus		+ KC_20_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TSE])  = (KC_20_0 * sqr_velz_plus_velZ	+ KC_20_1 * velz_term_plus		+ KC_20_2) * real_t(0.5);	 
	 src->get(x,y,z,Stencil_T::idx[E])    = KC_21_0 * oneminus_sqr_velz	- KC_21_1 * real_t(2.0) * velZ		- KC_21_2;
	 src->get(x,y,z,Stencil_T::idx[BE])   = (KC_21_0 * sqr_velz_minus_velZ	+ KC_21_1 * velz_term_minus		+ KC_21_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TE])   = (KC_21_0 * sqr_velz_plus_velZ	+ KC_21_1 * velz_term_plus		+ KC_21_2) * real_t(0.5);	 
	 src->get(x,y,z,Stencil_T::idx[NE])   = KC_22_0 * oneminus_sqr_velz	- KC_22_1 * real_t(2.0) * velZ		- KC_22_2;
	 src->get(x,y,z,Stencil_T::idx[BNE])  = (KC_22_0 * sqr_velz_minus_velZ	+ KC_22_1 * velz_term_minus		+ KC_22_2) * real_t(0.5);
	 src->get(x,y,z,Stencil_T::idx[TNE])  = (KC_22_0 * sqr_velz_plus_velZ	+ KC_22_1 * velz_term_plus		+ KC_22_2) * real_t(0.5);
         
	 WALBERLA_CHECK( math::finite( src->get(x,y,z,Stencil_T::idx[SW]) ) );
	 
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

// undefining the macros 
#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_CUMULANT_1

} // namespace lbm
} // namespace walberla 
