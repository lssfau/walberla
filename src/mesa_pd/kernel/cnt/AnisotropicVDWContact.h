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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <core/math/Angles.h>
#include <core/math/Constants.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
namespace cnt {

/**
 * anisotropic vdW contact
 */
class AnisotropicVDWContact
{
public:
   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor &ac);

   static constexpr real_t eps_ = 0.03733_r;
   static constexpr real_t A_ = 0.0223_r;
   static constexpr real_t B_ = 1.31_r;
   static constexpr real_t alf_ = 9.5_r;
   static constexpr real_t bet_ = 4.0_r;
   static constexpr real_t Cg_ = 90_r;
   static constexpr real_t del_ = -7.5_r;

   static constexpr size_t M = 5;
   /// Set of fitted constants, that determine the level of "smoothness" of vdW potential -
   /// magnitude of shear forces between sliding CNTs
   static constexpr std::array<real_t, M> C = { 0.35819_r, 0.03263_r, -0.00138_r, -0.00017_r, 0.00024_r };

   /// vdW interaction radius w/r to inertial segment radius
   static constexpr real_t CutoffFactor     = 4_r;
   /// CNT radius
   static constexpr real_t R_CNT = 6.78_r;
   static constexpr real_t R_ = R_CNT;
   static constexpr real_t r_cut_ = CutoffFactor * R_;

   auto isParallel() const {return isParallel_;}
   auto getLastEnergy() const {return energy_;}
private:
   real_t energy_; ///< total potential
   bool isParallel_;
};

template<typename Accessor>
inline
void AnisotropicVDWContact::operator()(const size_t p_idx1,
                                       const size_t p_idx2,
                                       Accessor &ac)
{
   //===Adaptation of PFC5 vdW contact model implementation====

   // Getting the orientations of segments
   Vec3 b1 = ac.getRotation(p_idx1).getMatrix() * Vec3(1.0, 0.0, 0.0); ///< ball 1 axial direction
   Vec3 b2 = ac.getRotation(p_idx2).getMatrix() * Vec3(1.0, 0.0, 0.0); ///< ball 2 axial direction


   // Distance between segments

   Vec3 n = ac.getPosition(p_idx2) - ac.getPosition(p_idx1); ///< contact normal
   auto L = n.length();
   n *= (1_r/L);


   //WALBERLA_LOG_DEVEL( "Normal: n = " << n );
   //WALBERLA_LOG_DEVEL( "Orientation of seg 2: b2 = " << b2 );
   //WALBERLA_LOG_DEVEL( "Orientation of seg 1: b1 = " << b1 );
   //WALBERLA_LOG_DEVEL( "Length of rad vect: L = " << L );


   constexpr real_t TOL = 10e-8_r;
   //---------------------
   // NORMALS CALCULATION
   //---------------------
   // names convention:
   // c1 - contact 1-2 normal
   // b1 - ball 1 axial direction
   // b2 - ball 2 axial direction
   // b3 - neytral direction
   // g - alighning torque direction
   // d - neytral plane normal direction
   // s - shear force direction

   // angle gamma - angle between two axial directions
   auto cos_gamma = b1 * b2;



   // if the angle between two axal directions is blunt, then inverce b2
   if (cos_gamma < 0_r)
   {
      b2 = -b2;
      cos_gamma = -cos_gamma;
   }
   // check that cosine belongs [-1,1]
   cos_gamma = std::min(1.0_r, cos_gamma);
   cos_gamma = std::max(-1.0_r, cos_gamma);
   isParallel_ = false;
   if (L < 20_r && L > 16_r)
   {
      const auto gamma = acos(cos_gamma);
      if (gamma < math::degToRad(10_r) || gamma > math::degToRad(170_r))
         isParallel_ = true;
   }
   //WALBERLA_LOG_DEVEL( "cos_gamma: = " << cos_gamma );


   // calculate functions of double argument
   auto sin_gamma = std::sqrt(1.0_r - cos_gamma * cos_gamma);
   auto cos_2gamma = cos_gamma * cos_gamma - sin_gamma * sin_gamma;
   auto sin_2gamma = 2.0_r * sin_gamma * cos_gamma;

   //WALBERLA_LOG_DEVEL( "sin_gamma: = " << sin_gamma );
   //WALBERLA_LOG_DEVEL( "cos_2gamma: = " << cos_2gamma );
   //WALBERLA_LOG_DEVEL( "sin_2gamma: = " << sin_2gamma );


   // g - direction of the aligning torques - b1 X b2
   Vec3 g(0.0, 0.0, 0.0);
   if (sin_gamma > TOL)
   {
      g = b1 % b2;
      g = g * (1.0_r / g.length());
   }
   //WALBERLA_LOG_DEVEL( "Aligning moment direction: g = " << g );

   // b3 - vector defining the neutral plane ( plane of shear forces )
   Vec3 b3 = b1 + b2;
   b3 = b3 * (1.0_r / b3.length());
   //WALBERLA_LOG_DEVEL( "Neutral plane defined by b3 = " << b3 );

   // angle theta - angle between b3 and c1
   auto cos_theta = b3 * n;
   // check that cosine belongs [-1,1]
   cos_theta = std::min(1.0_r, cos_theta);
   cos_theta = std::max(-1.0_r, cos_theta);
   //WALBERLA_LOG_DEVEL( "cos_theta: = " << cos_theta );

   // calculation of shear force direction
   Vec3 s(0.0, 0.0, 0.0);
   Vec3 d(0.0, 0.0, 0.0);

   if ((cos_theta > -1.0 + TOL) || (cos_theta < 1.0 - TOL))
      d = n % b3;
   s = n % d;
   s = s * (1.0_r / s.length());

   //WALBERLA_LOG_DEVEL( "Shear force direction: = " << s );
   //--------------------------------
   // NORMALS CALCULATION - END
   //--------------------------------

   // Fast calculation of trigonometric functions ( Chebyshev formulas )
   real_t coss[M], sinn[M];
   real_t sin_theta = std::sqrt(1.0_r - cos_theta * cos_theta);

   coss[0] = cos_theta * cos_theta - sin_theta * sin_theta;
   sinn[0] = 2.0_r * sin_theta * cos_theta;

   for (size_t i = 0; i < M - 1; ++i)
   {
      coss[i + 1] = coss[i] * coss[0] - sinn[i] * sinn[0];
      sinn[i + 1] = sinn[i] * coss[0] + sinn[0] * coss[i];
   }

   //WALBERLA_LOG_DEVEL( "coss: = " << coss[0] <<" "<< coss[1] <<" "<< coss[2] <<" "<< coss[3] <<" "<< coss[4]);
   //WALBERLA_LOG_DEVEL( "sinn: = " << sinn[0] <<" "<< sinn[1] <<" "<< sinn[2] <<" "<< sinn[3] <<" "<< sinn[4]);

   //WALBERLA_LOG_DEVEL( "C: = " << C[0] <<" "<< C[1] <<" "<< C[2] <<" "<< C[3] <<" "<< C[4]);

   // Cutoff for theta adjustment
   real_t W_th = 1_r;
   real_t W_th_L = 0_r;
   real_t W_th_LL = 0_r;

   // Adjustment w/r to O
   real_t TH = 1_r;
   real_t TH_L = 0_r;
   real_t TH_LL = 0_r;
   real_t TH_O = 0_r;
   real_t TH_OO = 0_r;

   real_t sign = 1_r;
   for (size_t i = 0; i < M; ++i)
   {
      TH = TH + W_th * C[i] * (sign + coss[i]);
      TH_L = TH_L + W_th_L * C[i] * (sign + coss[i]);
      TH_LL = TH_LL + W_th_LL * C[i] * (sign + coss[i]);
      TH_O = TH_O - W_th * 2_r * real_t(i + 1) * C[i] * (sinn[i]);
      TH_OO = TH_OO - W_th * 4_r * real_t(i + 1) * real_t(i + 1) * C[i] * (coss[i]);
      sign *= -1_r;
   }

   //WALBERLA_LOG_DEVEL( "TH: = " << TH );
   //WALBERLA_LOG_DEVEL( "TH_L: = " << TH_L );
   //WALBERLA_LOG_DEVEL( "TH_LL: = " << TH_LL );
   //WALBERLA_LOG_DEVEL( "TH_O: = " << TH_O );
   //WALBERLA_LOG_DEVEL( "TH_OO: = " << TH_OO );


   //------------------------------------------------------------------
   // THIS BLOCK IMPLEMENTS IF THE DISTANCE L IS WITHIN WORKING RANGE
   //------------------------------------------------------------------
   if ((L < 2_r * r_cut_) && (L > 2_r * R_ * 1.2_r * TH))
   {
      //-----Constants that appear in the potential--------------------------
      // This set of constants is described in our paper.
      //---------------------------------------------------------------------

      //-----Function D and its derivatives-----------------------
      real_t D = L / (R_ * TH) - 2_r;
      real_t D_L = 1_r / (R_ * TH) - (L * TH_L) / (R_ * TH * TH);
      real_t D_O = -(L * TH_O) / (R_ * TH * TH);
      real_t D_LL = -(TH_L) / (R_ * TH * TH);
      D_LL = D_LL - ((R_ * TH * TH) * (TH_L + L * TH_LL)) / (R_ * R_ * TH * TH * TH * TH);
      D_LL = D_LL + (2_r * R_ * L * TH * TH_L * TH_L) / (R_ * R_ * TH * TH * TH * TH);
      real_t D_OO = -(R_ * L * TH_OO * TH * TH) / (R_ * R_ * TH * TH * TH * TH);
      D_OO = D_OO + (2_r * R_ * L * TH * TH_O * TH_O) / (R_ * R_ * TH * TH * TH * TH);
      //-----------------------------------------------------------

      //WALBERLA_LOG_DEVEL( "D: = " << D );
      //WALBERLA_LOG_DEVEL( "D_L: = " << D_L );
      //WALBERLA_LOG_DEVEL( "D_LL: = " << D_LL );
      //WALBERLA_LOG_DEVEL( "D_O: = " << D_O );
      //WALBERLA_LOG_DEVEL( "D_OO: = " << D_OO );


      //----Function Vc and its derivatives---------------------------------------
      const real_t DpowAlpha0 = std::pow(D, -(alf_));
      const real_t DpowAlpha1 = std::pow(D, -(alf_ + 1));
      const real_t DpowAlpha2 = std::pow(D, -(alf_ + 2));
      const real_t DpowBeta0 = std::pow(D, -(bet_));
      const real_t DpowBeta1 = std::pow(D, -(bet_ + 1));
      const real_t DpowBeta2 = std::pow(D, -(bet_ + 2));
      real_t Vc = 4_r * eps_ * (A_ * DpowAlpha0 - B_ * DpowBeta0);
      real_t Vc_L = 4_r * eps_ * (-alf_ * A_ * DpowAlpha1 + bet_ * B_ * DpowBeta1) * D_L;
      real_t Vc_O = 4_r * eps_ * (-alf_ * A_ * DpowAlpha1 + bet_ * B_ * DpowBeta1) * D_O;
      real_t Vc_LL = 4_r * eps_ * (alf_ * (alf_ + 1_r) * A_ * DpowAlpha2 - bet_ * (bet_ + 1_r) * B_ * DpowBeta2) * D_L;
      Vc_LL = Vc_LL + 4_r * eps_ * (-alf_ * A_ * DpowAlpha1 + bet_ * B_ * DpowBeta1) * D_LL;
      real_t Vc_OO = 4_r * eps_ * (alf_ * (alf_ + 1_r) * A_ * DpowAlpha2 - bet_ * (bet_ + 1_r) * B_ * DpowBeta2) * D_O;
      Vc_OO = Vc_OO + 4_r * eps_ * (-alf_ * A_ * DpowAlpha1 + bet_ * B_ * DpowBeta1) * D_OO;
      //--------------------------------------------------------------------------

      //WALBERLA_LOG_DEVEL( "VC = " << Vc );
      //WALBERLA_LOG_DEVEL( "VC_L = " << Vc_L );
      //WALBERLA_LOG_DEVEL( "VC_LL = " << Vc_LL );

      //WALBERLA_LOG_DEVEL( "VC_O = " << Vc_O );
      //WALBERLA_LOG_DEVEL( "VC_OO = " << Vc_OO );


      // Cutoff for u adjustment
      real_t W_u = 1_r;
      real_t W_u_L = 0_r;
      real_t W_u_LL = 0_r;
      WALBERLA_UNUSED(W_u_LL);

      // Cubic cutoff function 3T->4T (hardcoded since we do not need to mess w these parameters)
      constexpr auto Q1_ = -80.0_r;
      constexpr auto Q2_ = 288.0_r;
      constexpr auto Q3_ = -336.0_r;
      constexpr auto Q4_ = 128.0_r;
      const real_t rcut2inv = 1_r / (2.0_r * r_cut_);
      real_t nd = L * rcut2inv;
      if ((nd > 0.75_r) && (nd < 1.0_r))
      {
         W_u = Q1_ + Q2_ * nd + Q3_ * nd * nd + Q4_ * nd * nd * nd;
         W_u_L = (Q2_ + 2.0_r * Q3_ * nd + 3.0_r * Q4_ * nd * nd) * rcut2inv;
         W_u_LL = (2.0_r * Q3_ + 6.0_r * Q4_ * nd) * rcut2inv * rcut2inv;
      }
      //--------------------------------------------------------------------------

      //WALBERLA_LOG_DEVEL( "W_u = " << W_u );
      //WALBERLA_LOG_DEVEL( "W_u_L = " << W_u_L );
      //WALBERLA_LOG_DEVEL( "W_u_LL = " << W_u_LL );


      // Cutoff for gamma adjustment

      real_t W_ga, W_ga_L, W_ga_LL;
      if (L / R_ > 2.75_r)
      {
         W_ga = Cg_ * std::pow((L / R_), del_);
         W_ga_L = ((del_ * Cg_) / R_) * std::pow((L / R_), del_ - 1_r);
         W_ga_LL = ((del_ * (del_ - 1) * Cg_) / (R_ * R_)) * std::pow((L / R_), del_ - 2_r);
      } else
      {
         W_ga = Cg_ * std::pow((2.75_r), del_);
         W_ga_L = 0;
         W_ga_LL = 0;
      }

      //WALBERLA_LOG_DEVEL( "W_ga = " << W_ga );
      //WALBERLA_LOG_DEVEL( "W_ga_L = " << W_ga_L );
      //WALBERLA_LOG_DEVEL( "W_ga_LL = " << W_ga_LL );


      real_t GA = 1_r;
      real_t GA_L = 0_r;
      real_t GA_LL = 0_r;
      WALBERLA_UNUSED(GA_LL);
      real_t GA_G = 0_r;
      real_t GA_GG = 0_r;
      WALBERLA_UNUSED(GA_GG);
      if (std::abs(sin_gamma) > TOL)
      {
         GA = 1_r + W_ga * (1_r - cos_2gamma);
         GA_L = W_ga_L * (1_r - cos_2gamma);
         GA_LL = W_ga_LL * (1_r - cos_2gamma);
         GA_G = 2_r * W_ga * sin_2gamma;
         GA_GG = 4_r * W_ga * cos_2gamma;
      }

      //----Forces and torque-----------------------
      real_t FL = -GA_L * W_u * Vc - GA * W_u_L * Vc - GA * W_u * Vc_L;
      real_t FO = -(1_r / L) * GA * W_u * Vc_O;
      real_t MG = -GA_G * W_u * Vc;

      //WALBERLA_LOG_DEVEL( "FL = " << FL );
      //WALBERLA_LOG_DEVEL( "FO = " << FO );
      //WALBERLA_LOG_DEVEL( "MG = " << MG );


      Vec3 force = FL * n + FO * s;
      Vec3 moment = MG * g;


      //WALBERLA_LOG_DEVEL( "Contact force: = " << force );
      //WALBERLA_LOG_DEVEL( "Contact moment: = " << moment );

      addForceAtomic(p_idx1, ac, -force);
      addForceAtomic(p_idx2, ac,  force);

      addTorqueAtomic(p_idx1, ac, -moment);
      addTorqueAtomic(p_idx2, ac,  moment);

      // Potential energy
      energy_ = GA * W_u * Vc;
      // WALBERLA_LOG_DEVEL( "U_vdw = " << U );
   } else if (L <= 2_r * R_ * 1.2_r * TH)
   { // Small distance
      //WALBERLA_LOG_DEVEL( "Small distance");
      energy_ = 0_r;
      real_t F = -1_r;
      addForceAtomic(p_idx1, ac,  F * n);
      addForceAtomic(p_idx2, ac, -F * n);
   }
}

} //namespace cnt
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla