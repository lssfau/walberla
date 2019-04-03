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
//! \author Ehsan Fattahi <ehsan.fattahi@fau.de>
//! \author Felix Winterhalter <felix.winterhalter@fau.de>
//
//======================================================================================================================

#include "lbm/field/DensityVelocityCallback.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"
#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include <type_traits>


namespace walberla {
namespace lbm {

//////////////////////////
// D3Q19 SPECIALIZATION //
//////////////////////////

// This MRT variant is taken from the dissertation of Ulf Schiller 
// "Thermal fluctuations and boundary conditions in the lattice Boltzmann method" (2008), p. 24ff
// There are some typos in the moment matrix on p.27
// The here implemented ordering of the moments is however different from the reference (Eq. 2.61-2.63)
// The moments are weighted orthogonal (Eq. 2.58)

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 \
   std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >::value && \
   std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value && \
   ! LatticeModel_T::compressible && \
   LatticeModel_T::equilibriumAccuracyOrder == 2 && \
   std::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation >::value

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )
{
   const auto & collisionModel = src->latticeModel().collisionModel();

   const real_t s1  = collisionModel.s1();
   const real_t s2  = collisionModel.s2();
   const real_t s4  = collisionModel.s4();
   const real_t s6  = collisionModel.s6();
   const real_t s8  = collisionModel.s8();
   const real_t s9  = collisionModel.s9();
   const real_t s10 = collisionModel.s10();
   const real_t s11 = collisionModel.s11();
   const real_t s12 = collisionModel.s12();
   const real_t s13 = collisionModel.s13();
   const real_t s14 = collisionModel.s14();
   const real_t s15 = collisionModel.s15();
   const real_t s16 = collisionModel.s16();
   const real_t s17 = collisionModel.s17();
   const real_t s18 = collisionModel.s18();

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()

         const Vector3<real_t> velocity( velX, velY, velZ );
         this->densityVelocityOut( x, y, z, lm, velocity, rho + real_t(1) );

         const real_t velSqr = velX * velX + velY * velY + velZ * velZ;

         const real_t vel9 = real_t(2) * velX * velX - velY * velY - velZ * velZ;
         const real_t vel11 = velY * velY - velZ * velZ;
         const real_t vel13 = velX * velY;
         const real_t vel14 = velY * velZ;
         const real_t vel15 = velX * velZ;

         const real_t mStar0  = rho;
         const real_t mStar1  = velSqr + ( real_t(1) - s1 ) * ( -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE - velSqr );
         const real_t mStar2  = ( real_t(1) - s2 ) * ( vC - real_t(2) * ( vN + vS + vW + vE + vT + vB )
                                                          + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE );
         const real_t mStar3  = velX;
         const real_t mStar4  = ( real_t(1) - s4 ) * ( real_t(2) * vW - real_t(2) * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE );
         const real_t mStar5  = velY;
         const real_t mStar6  = ( real_t(1) - s6 ) * ( real_t(-2) * vN + real_t(2) * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS );
         const real_t mStar7  = velZ;
         const real_t mStar8  = ( real_t(1) - s8 ) * ( real_t(-2) * vT + real_t(2) * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE );
         const real_t mStar9  = vel9 + ( real_t(1) - s9 ) * ( -vN - vS + real_t(2) * vW + real_t(2) * vE - vT - vB + vNW + vNE + vSW + vSE - real_t(2) * vTN -
                                                              real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE - vel9 );
         const real_t mStar10 = ( real_t(1) - s10 ) * ( vN + vS - real_t(2) * vW - real_t(2) * vE + vT + vB + vNW + vNE + vSW + vSE - real_t(2) * vTN -
                                                        real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE );
         const real_t mStar11 = vel11 + ( real_t(1) - s11 ) * ( vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE - vel11 );
         const real_t mStar12 = ( real_t(1) - s12 ) * ( -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE );
         const real_t mStar13 = vel13 + ( real_t(1) - s13 ) * ( -vNW + vNE + vSW - vSE - vel13 );
         const real_t mStar14 = vel14 + ( real_t(1) - s14 ) * (  vTN - vTS - vBN + vBS - vel14 );
         const real_t mStar15 = vel15 + ( real_t(1) - s15 ) * ( -vTW + vTE + vBW - vBE - vel15 );
         const real_t mStar16 = ( real_t(1) - s16 ) * ( -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE );
         const real_t mStar17 = ( real_t(1) - s17 ) * ( -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS );
         const real_t mStar18 = ( real_t(1) - s18 ) * ( -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE );

         dst->get( x, y, z, Stencil_T::idx[C] )  = _1_3  * mStar0  - _1_2  * mStar1  + _1_6  * mStar2;
         dst->get( x, y, z, Stencil_T::idx[N] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar5  - _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[S] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar5  + _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[W] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar3  + _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         dst->get( x, y, z, Stencil_T::idx[E] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar3  - _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         dst->get( x, y, z, Stencil_T::idx[T] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar7  - _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[B] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar7  + _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 - _1_8  * mStar16 - _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[NE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 + _1_8  * mStar16 - _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 - _1_8  * mStar16 + _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[SE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 + _1_8  * mStar16 + _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[TN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 +
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 -
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 + _1_8  * mStar16 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 - _1_8  * mStar16 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 +
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 -
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 + _1_8  * mStar16 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 - _1_8  * mStar16 - _1_8  * mStar18;
         
         if (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value == false)
         {
            const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho + real_t(1.0), collisionModel.omega(), collisionModel.omega_bulk() );
            for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
               dst->get( x, y, z, d.toIdx() ) += lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho + real_t(1.0), commonForceTerms, LatticeModel_T::w[ d.toIdx() ], real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), collisionModel.omega(), collisionModel.omega_bulk() );
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )
{
   const auto & collisionModel = src->latticeModel().collisionModel();

   const real_t s1  = collisionModel.s1();
   const real_t s2  = collisionModel.s2();
   const real_t s4  = collisionModel.s4();
   const real_t s6  = collisionModel.s6();
   const real_t s8  = collisionModel.s8();
   const real_t s9  = collisionModel.s9();
   const real_t s10 = collisionModel.s10();
   const real_t s11 = collisionModel.s11();
   const real_t s12 = collisionModel.s12();
   const real_t s13 = collisionModel.s13();
   const real_t s14 = collisionModel.s14();
   const real_t s15 = collisionModel.s15();
   const real_t s16 = collisionModel.s16();
   const real_t s17 = collisionModel.s17();
   const real_t s18 = collisionModel.s18();

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

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         const Vector3<real_t> velocity( velX, velY, velZ );
         this->densityVelocityOut( x, y, z, lm, velocity, rho + real_t(1) );

         const real_t velSqr = velX * velX + velY * velY + velZ * velZ;

         const real_t vel9 = real_t(2) * velX * velX - velY * velY - velZ * velZ;
         const real_t vel11 = velY * velY - velZ * velZ;
         const real_t vel13 = velX * velY;
         const real_t vel14 = velY * velZ;
         const real_t vel15 = velX * velZ;

         const real_t mStar0  = rho;
         const real_t mStar1  = velSqr + ( real_t(1) - s1 ) * ( -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE - velSqr );
         const real_t mStar2  = ( real_t(1) - s2 ) * ( vC - real_t(2) * ( vN + vS + vW + vE + vT + vB )
                                                          + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE );
         const real_t mStar3  = velX;
         const real_t mStar4  = ( real_t(1) - s4 ) * ( real_t(2) * vW - real_t(2) * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE );
         const real_t mStar5  = velY;
         const real_t mStar6  = ( real_t(1) - s6 ) * ( real_t(-2) * vN + real_t(2) * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS );
         const real_t mStar7  = velZ;
         const real_t mStar8  = ( real_t(1) - s8 ) * ( real_t(-2) * vT + real_t(2) * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE );
         const real_t mStar9  = vel9 + ( real_t(1) - s9 ) * ( -vN - vS + real_t(2) * vW + real_t(2) * vE - vT - vB + vNW + vNE + vSW + vSE - real_t(2) * vTN -
                                                              real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE - vel9 );
         const real_t mStar10 = ( real_t(1) - s10 ) * ( vN + vS - real_t(2) * vW - real_t(2) * vE + vT + vB + vNW + vNE + vSW + vSE - real_t(2) * vTN -
                                                        real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE );
         const real_t mStar11 = vel11 + ( real_t(1) - s11 ) * ( vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE - vel11 );
         const real_t mStar12 = ( real_t(1) - s12 ) * ( -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE );
         const real_t mStar13 = vel13 + ( real_t(1) - s13 ) * ( -vNW + vNE + vSW - vSE - vel13 );
         const real_t mStar14 = vel14 + ( real_t(1) - s14 ) * (  vTN - vTS - vBN + vBS - vel14 );
         const real_t mStar15 = vel15 + ( real_t(1) - s15 ) * ( -vTW + vTE + vBW - vBE - vel15 );
         const real_t mStar16 = ( real_t(1) - s16 ) * ( -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE );
         const real_t mStar17 = ( real_t(1) - s17 ) * ( -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS );
         const real_t mStar18 = ( real_t(1) - s18 ) * ( -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE );

         src->get( x, y, z, Stencil_T::idx[C] )  = _1_3  * mStar0  - _1_2  * mStar1  + _1_6  * mStar2;
         src->get( x, y, z, Stencil_T::idx[N] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar5  - _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[S] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar5  + _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[W] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar3  + _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         src->get( x, y, z, Stencil_T::idx[E] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar3  - _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         src->get( x, y, z, Stencil_T::idx[T] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar7  - _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[B] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar7  + _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[NW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 - _1_8  * mStar16 - _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[NE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 + _1_8  * mStar16 - _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[SW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 - _1_8  * mStar16 + _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[SE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 + _1_8  * mStar16 + _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[TN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 +
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 -
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 + _1_8  * mStar16 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 - _1_8  * mStar16 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 +
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 -
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 + _1_8  * mStar16 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 - _1_8  * mStar16 - _1_8  * mStar18;
         
         if (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value == false)
         {
            const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho + real_t(1.0), collisionModel.omega(), collisionModel.omega_bulk() );
            for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
               src->get( x, y, z, d.toIdx() ) += lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho + real_t(1.0), commonForceTerms, LatticeModel_T::w[ d.toIdx() ], real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), collisionModel.omega(), collisionModel.omega_bulk() );
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1



} // namespace lbm
} // namespace walberla
