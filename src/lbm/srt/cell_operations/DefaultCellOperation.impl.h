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
//! \file DefaultCellOperation.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"

#include <type_traits>


namespace walberla {
namespace lbm {



////////////////////////////////////////////////////////////////////////
// Available SRT implementations:                                     //
//                                                                    //
// Generic (D*Q*) version:                                            //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
//                                                                    //
// Optimized D3Q19 implementation:                                    //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
////////////////////////////////////////////////////////////////////////


///////////////////////////////
// Specialization for:       //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                     LatticeModel_T::CollisionModel::constant &&
                                                                     std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false),              "There is a specialization for D3Q19!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   DefaultCellOperation() : omega_( real_t(0) ), latticeModel_( NULL ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      omega_ = latticeModel.collisionModel().omega();
      latticeModel_ = &latticeModel;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = dst->getDensityAndVelocity( velocity, x, y, z );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
                                                          omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho );
      }
   }

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = getDensityAndVelocity( velocity, *latticeModel_, dst );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst[ d.toIdx() ] = ( real_t(1.0) - omega_ ) * dst[ d.toIdx() ] +
                                            omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho );
      }
   }
private:

   real_t omega_;
   const LatticeModel_T * latticeModel_;
};



/////////////////////////////////////////////
// Specialization for:                     //
// - compressible                          //
// - additional forces (Guo - constant)    //
/////////////////////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                     LatticeModel_T::CollisionModel::constant &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Guo_tag >::value
                                                                     >::type >
{                                                                                   
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Guo_tag >::value),         "Only works with Guo constant force model !" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   DefaultCellOperation() : omega_( real_t(0) ), latticeModel_( NULL ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      omega_ = latticeModel.collisionModel().omega();
      latticeModel_ = &latticeModel;
      force_ = latticeModel.forceModel().forceDensity();
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = dst->getDensityAndVelocity( velocity, x, y, z );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );

         const real_t force_trm = real_t(3.0) * LatticeModel_T::w[ d.toIdx() ] * ( real_t(1) - real_t(0.5) * omega_ ) *
                                  ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * force_ );

         dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
                                                          omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                                          force_trm;
      }
   }

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = getDensityAndVelocity( velocity, *latticeModel_, dst );

      // collide< LatticeModel_T::CollisionModel::constant >
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );

         const real_t force_trm = real_t(3.0) * LatticeModel_T::w[ d.toIdx() ] * ( real_t(1) - real_t(0.5) * omega_ ) *
                                  ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * force_ );

         dst[ d.toIdx() ] = ( real_t(1.0) - omega_ ) * dst[ d.toIdx() ] +
                                            omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                              force_trm;
      }
   }

private:

   real_t omega_;
   const LatticeModel_T * latticeModel_;
   Vector3<real_t> force_;
};







///////////////////////////////
// Specialization for D3Q19: //
// - incompressible          //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                     LatticeModel_T::CollisionModel::constant &&
                                                                     std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     ! LatticeModel_T::compressible &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   DefaultCellOperation() : omega_trm_( real_t(0) ), omega_w0_( real_t(0) ), omega_w1_( real_t(0) ), omega_w2_( real_t(0) ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      const real_t omega = latticeModel.collisionModel().omega();
      omega_trm_ = real_t(1) - omega;
      omega_w0_  = real_t(3) * ( real_t(1) / real_t( 3) ) * omega;
      omega_w1_  = real_t(3) * ( real_t(1) / real_t(18) ) * omega;
      omega_w2_  = real_t(3) * ( real_t(1) / real_t(36) ) * omega;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;

private:

   real_t omega_trm_;
   real_t omega_w0_;
   real_t omega_w1_;
   real_t omega_w2_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                    LatticeModel_T::CollisionModel::constant &&
                                                                    std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    ! LatticeModel_T::compressible &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   using namespace stencil;

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

   const real_t velX = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
   const real_t velY = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
   const real_t velZ = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) * rho - real_t(0.5) * ( velXX + velYY + velZZ );

   dst->get(x,y,z,Stencil::idx[C]) = omega_trm_ * dd_tmp_C + omega_w0_ * dir_indep_trm;

   const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

   dst->get(x,y,z,Stencil::idx[E]) = omega_trm_ * dd_tmp_E + omega_w1_ * ( vel_trm_E_W + velX );
   dst->get(x,y,z,Stencil::idx[W]) = omega_trm_ * dd_tmp_W + omega_w1_ * ( vel_trm_E_W - velX );
   dst->get(x,y,z,Stencil::idx[N]) = omega_trm_ * dd_tmp_N + omega_w1_ * ( vel_trm_N_S + velY );
   dst->get(x,y,z,Stencil::idx[S]) = omega_trm_ * dd_tmp_S + omega_w1_ * ( vel_trm_N_S - velY );
   dst->get(x,y,z,Stencil::idx[T]) = omega_trm_ * dd_tmp_T + omega_w1_ * ( vel_trm_T_B + velZ );
   dst->get(x,y,z,Stencil::idx[B]) = omega_trm_ * dd_tmp_B + omega_w1_ * ( vel_trm_T_B - velZ );

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

   dst->get(x,y,z,Stencil::idx[NW]) = omega_trm_ * dd_tmp_NW + omega_w2_ * ( vel_trm_NW_SE - velXmY );
   dst->get(x,y,z,Stencil::idx[SE]) = omega_trm_ * dd_tmp_SE + omega_w2_ * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

   dst->get(x,y,z,Stencil::idx[NE]) = omega_trm_ * dd_tmp_NE + omega_w2_ * ( vel_trm_NE_SW + velXpY );
   dst->get(x,y,z,Stencil::idx[SW]) = omega_trm_ * dd_tmp_SW + omega_w2_ * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

   dst->get(x,y,z,Stencil::idx[TW]) = omega_trm_ * dd_tmp_TW + omega_w2_ * ( vel_trm_TW_BE - velXmZ );
   dst->get(x,y,z,Stencil::idx[BE]) = omega_trm_ * dd_tmp_BE + omega_w2_ * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

   dst->get(x,y,z,Stencil::idx[TE]) = omega_trm_ * dd_tmp_TE + omega_w2_ * ( vel_trm_TE_BW + velXpZ );
   dst->get(x,y,z,Stencil::idx[BW]) = omega_trm_ * dd_tmp_BW + omega_w2_ * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

   dst->get(x,y,z,Stencil::idx[TS]) = omega_trm_ * dd_tmp_TS + omega_w2_ * ( vel_trm_TS_BN - velYmZ );
   dst->get(x,y,z,Stencil::idx[BN]) = omega_trm_ * dd_tmp_BN + omega_w2_ * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

   dst->get(x,y,z,Stencil::idx[TN]) = omega_trm_ * dd_tmp_TN + omega_w2_ * ( vel_trm_TN_BS + velYpZ );
   dst->get(x,y,z,Stencil::idx[BS]) = omega_trm_ * dd_tmp_BS + omega_w2_ * ( vel_trm_TN_BS - velYpZ );
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                    LatticeModel_T::CollisionModel::constant &&
                                                                    std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    ! LatticeModel_T::compressible &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src.neighbor(-1, -1,  0, Stencil::idx[NE]);
   const real_t dd_tmp_N  = src.neighbor( 0, -1,  0, Stencil::idx[N] );
   const real_t dd_tmp_NW = src.neighbor(+1, -1,  0, Stencil::idx[NW]);
   const real_t dd_tmp_W  = src.neighbor(+1,  0,  0, Stencil::idx[W] );
   const real_t dd_tmp_SW = src.neighbor(+1, +1,  0, Stencil::idx[SW]);
   const real_t dd_tmp_S  = src.neighbor( 0, +1,  0, Stencil::idx[S] );
   const real_t dd_tmp_SE = src.neighbor(-1, +1,  0, Stencil::idx[SE]);
   const real_t dd_tmp_E  = src.neighbor(-1,  0,  0, Stencil::idx[E] );
   const real_t dd_tmp_T  = src.neighbor( 0,  0, -1, Stencil::idx[T] );
   const real_t dd_tmp_TE = src.neighbor(-1,  0, -1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src.neighbor( 0, -1, -1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src.neighbor(+1,  0, -1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src.neighbor( 0, +1, -1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src.neighbor( 0,  0, +1, Stencil::idx[B] );
   const real_t dd_tmp_BE = src.neighbor(-1,  0, +1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src.neighbor( 0, -1, +1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src.neighbor(+1,  0, +1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src.neighbor( 0, +1, +1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src.neighbor( 0,  0,  0, Stencil::idx[C] );

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

   const real_t velX = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
   const real_t velY = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
   const real_t velZ = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) * rho - real_t(0.5) * ( velXX + velYY + velZZ );

   dst[ Stencil::idx[C] ] = omega_trm_ * dd_tmp_C + omega_w0_ * dir_indep_trm;

   const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

   dst[ Stencil::idx[E] ] = omega_trm_ * dd_tmp_E + omega_w1_ * ( vel_trm_E_W + velX );
   dst[ Stencil::idx[W] ] = omega_trm_ * dd_tmp_W + omega_w1_ * ( vel_trm_E_W - velX );
   dst[ Stencil::idx[N] ] = omega_trm_ * dd_tmp_N + omega_w1_ * ( vel_trm_N_S + velY );
   dst[ Stencil::idx[S] ] = omega_trm_ * dd_tmp_S + omega_w1_ * ( vel_trm_N_S - velY );
   dst[ Stencil::idx[T] ] = omega_trm_ * dd_tmp_T + omega_w1_ * ( vel_trm_T_B + velZ );
   dst[ Stencil::idx[B] ] = omega_trm_ * dd_tmp_B + omega_w1_ * ( vel_trm_T_B - velZ );

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

   dst[ Stencil::idx[NW] ] = omega_trm_ * dd_tmp_NW + omega_w2_ * ( vel_trm_NW_SE - velXmY );
   dst[ Stencil::idx[SE] ] = omega_trm_ * dd_tmp_SE + omega_w2_ * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

   dst[ Stencil::idx[NE] ] = omega_trm_ * dd_tmp_NE + omega_w2_ * ( vel_trm_NE_SW + velXpY );
   dst[ Stencil::idx[SW] ] = omega_trm_ * dd_tmp_SW + omega_w2_ * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

   dst[ Stencil::idx[TW] ] = omega_trm_ * dd_tmp_TW + omega_w2_ * ( vel_trm_TW_BE - velXmZ );
   dst[ Stencil::idx[BE] ] = omega_trm_ * dd_tmp_BE + omega_w2_ * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

   dst[ Stencil::idx[TE] ] = omega_trm_ * dd_tmp_TE + omega_w2_ * ( vel_trm_TE_BW + velXpZ );
   dst[ Stencil::idx[BW] ] = omega_trm_ * dd_tmp_BW + omega_w2_ * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

   dst[ Stencil::idx[TS] ] = omega_trm_ * dd_tmp_TS + omega_w2_ * ( vel_trm_TS_BN - velYmZ );
   dst[ Stencil::idx[BN] ] = omega_trm_ * dd_tmp_BN + omega_w2_ * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

   dst[ Stencil::idx[TN] ] = omega_trm_ * dd_tmp_TN + omega_w2_ * ( vel_trm_TN_BS + velYpZ );
   dst[ Stencil::idx[BS] ] = omega_trm_ * dd_tmp_BS + omega_w2_ * ( vel_trm_TN_BS - velYpZ );
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                     LatticeModel_T::CollisionModel::constant &&
                                                                     std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     LatticeModel_T::compressible &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   DefaultCellOperation() : omega_trm_( real_t(0) ), omega_w0_( real_t(0) ), omega_w1_( real_t(0) ), omega_w2_( real_t(0) ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      const real_t omega = latticeModel.collisionModel().omega();
      omega_trm_ = real_t(1) - omega;
      omega_w0_  = real_t(3) * ( real_t(1) / real_t( 3) ) * omega;
      omega_w1_  = real_t(3) * ( real_t(1) / real_t(18) ) * omega;
      omega_w2_  = real_t(3) * ( real_t(1) / real_t(36) ) * omega;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;


private:

   real_t omega_trm_;
   real_t omega_w0_;
   real_t omega_w1_;
   real_t omega_w2_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                    LatticeModel_T::CollisionModel::constant &&
                                                                    std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    LatticeModel_T::compressible &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   using namespace stencil;

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

   const real_t velX = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
   const real_t velY = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
   const real_t velZ = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) - real_t(0.5) * ( velXX + velYY + velZZ );

   dst->get(x,y,z,Stencil::idx[C]) = omega_trm_ * dd_tmp_C + omega_w0_ * rho * dir_indep_trm;

   const real_t omega_w1_rho = omega_w1_ * rho;

   const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

   dst->get(x,y,z,Stencil::idx[E]) = omega_trm_ * dd_tmp_E + omega_w1_rho * ( vel_trm_E_W + velX );
   dst->get(x,y,z,Stencil::idx[W]) = omega_trm_ * dd_tmp_W + omega_w1_rho * ( vel_trm_E_W - velX );
   dst->get(x,y,z,Stencil::idx[N]) = omega_trm_ * dd_tmp_N + omega_w1_rho * ( vel_trm_N_S + velY );
   dst->get(x,y,z,Stencil::idx[S]) = omega_trm_ * dd_tmp_S + omega_w1_rho * ( vel_trm_N_S - velY );
   dst->get(x,y,z,Stencil::idx[T]) = omega_trm_ * dd_tmp_T + omega_w1_rho * ( vel_trm_T_B + velZ );
   dst->get(x,y,z,Stencil::idx[B]) = omega_trm_ * dd_tmp_B + omega_w1_rho * ( vel_trm_T_B - velZ );

   const real_t omega_w2_rho = omega_w2_ * rho;

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

   dst->get(x,y,z,Stencil::idx[NW]) = omega_trm_ * dd_tmp_NW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
   dst->get(x,y,z,Stencil::idx[SE]) = omega_trm_ * dd_tmp_SE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

   dst->get(x,y,z,Stencil::idx[NE]) = omega_trm_ * dd_tmp_NE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
   dst->get(x,y,z,Stencil::idx[SW]) = omega_trm_ * dd_tmp_SW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

   dst->get(x,y,z,Stencil::idx[TW]) = omega_trm_ * dd_tmp_TW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
   dst->get(x,y,z,Stencil::idx[BE]) = omega_trm_ * dd_tmp_BE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

   dst->get(x,y,z,Stencil::idx[TE]) = omega_trm_ * dd_tmp_TE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
   dst->get(x,y,z,Stencil::idx[BW]) = omega_trm_ * dd_tmp_BW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

   dst->get(x,y,z,Stencil::idx[TS]) = omega_trm_ * dd_tmp_TS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
   dst->get(x,y,z,Stencil::idx[BN]) = omega_trm_ * dd_tmp_BN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

   dst->get(x,y,z,Stencil::idx[TN]) = omega_trm_ * dd_tmp_TN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
   dst->get(x,y,z,Stencil::idx[BS]) = omega_trm_ * dd_tmp_BS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                    LatticeModel_T::CollisionModel::constant &&
                                                                    std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    LatticeModel_T::compressible &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src.neighbor(-1, -1,  0, Stencil::idx[NE]);
   const real_t dd_tmp_N  = src.neighbor( 0, -1,  0, Stencil::idx[N] );
   const real_t dd_tmp_NW = src.neighbor(+1, -1,  0, Stencil::idx[NW]);
   const real_t dd_tmp_W  = src.neighbor(+1,  0,  0, Stencil::idx[W] );
   const real_t dd_tmp_SW = src.neighbor(+1, +1,  0, Stencil::idx[SW]);
   const real_t dd_tmp_S  = src.neighbor( 0, +1,  0, Stencil::idx[S] );
   const real_t dd_tmp_SE = src.neighbor(-1, +1,  0, Stencil::idx[SE]);
   const real_t dd_tmp_E  = src.neighbor(-1,  0,  0, Stencil::idx[E] );
   const real_t dd_tmp_T  = src.neighbor( 0,  0, -1, Stencil::idx[T] );
   const real_t dd_tmp_TE = src.neighbor(-1,  0, -1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src.neighbor( 0, -1, -1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src.neighbor(+1,  0, -1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src.neighbor( 0, +1, -1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src.neighbor( 0,  0, +1, Stencil::idx[B] );
   const real_t dd_tmp_BE = src.neighbor(-1,  0, +1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src.neighbor( 0, -1, +1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src.neighbor(+1,  0, +1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src.neighbor( 0, +1, +1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src.neighbor( 0,  0,  0, Stencil::idx[C] );

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;
   const real_t invRho = real_t(1.0) / rho;

   const real_t velX = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
   const real_t velY = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
   const real_t velZ = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) - real_t(0.5) * ( velXX + velYY + velZZ );

   dst[ Stencil::idx[C] ] = omega_trm_ * dd_tmp_C + omega_w0_ * rho * dir_indep_trm;

   const real_t omega_w1_rho = omega_w1_ * rho;

   const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

   dst[ Stencil::idx[E] ] = omega_trm_ * dd_tmp_E + omega_w1_rho * ( vel_trm_E_W + velX );
   dst[ Stencil::idx[W] ] = omega_trm_ * dd_tmp_W + omega_w1_rho * ( vel_trm_E_W - velX );
   dst[ Stencil::idx[N] ] = omega_trm_ * dd_tmp_N + omega_w1_rho * ( vel_trm_N_S + velY );
   dst[ Stencil::idx[S] ] = omega_trm_ * dd_tmp_S + omega_w1_rho * ( vel_trm_N_S - velY );
   dst[ Stencil::idx[T] ] = omega_trm_ * dd_tmp_T + omega_w1_rho * ( vel_trm_T_B + velZ );
   dst[ Stencil::idx[B] ] = omega_trm_ * dd_tmp_B + omega_w1_rho * ( vel_trm_T_B - velZ );

   const real_t omega_w2_rho = omega_w2_ * rho;

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

   dst[ Stencil::idx[NW] ] = omega_trm_ * dd_tmp_NW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
   dst[ Stencil::idx[SE] ] = omega_trm_ * dd_tmp_SE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

   dst[ Stencil::idx[NE] ] = omega_trm_ * dd_tmp_NE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
   dst[ Stencil::idx[SW] ] = omega_trm_ * dd_tmp_SW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

   dst[ Stencil::idx[TW] ] = omega_trm_ * dd_tmp_TW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
   dst[ Stencil::idx[BE] ] = omega_trm_ * dd_tmp_BE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

   dst[ Stencil::idx[TE] ] = omega_trm_ * dd_tmp_TE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
   dst[ Stencil::idx[BW] ] = omega_trm_ * dd_tmp_BW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

   dst[ Stencil::idx[TS] ] = omega_trm_ * dd_tmp_TS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
   dst[ Stencil::idx[BN] ] = omega_trm_ * dd_tmp_BN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

   dst[ Stencil::idx[TN] ] = omega_trm_ * dd_tmp_TN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
   dst[ Stencil::idx[BS] ] = omega_trm_ * dd_tmp_BS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
}












} // namespace lbm
} // namespace walberla
