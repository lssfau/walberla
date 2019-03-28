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
// Available TRT implementations:                                     //
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
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                     ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value), "Only works with TRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false),              "There is a specialization for D3Q19!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() : lambda_e_( real_t(0) ), lambda_d_( real_t(0) ), latticeModel_( NULL ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      lambda_e_ = latticeModel.collisionModel().lambda_e();
      lambda_d_ = latticeModel.collisionModel().lambda_d();
      latticeModel_ = &latticeModel;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;

private:

   real_t lambda_e_;
   real_t lambda_d_;
   const LatticeModel_T * latticeModel_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                    ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   real_t pdfs[ Stencil::Size ];
   
   // stream pull
   for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
   {
      const auto pdf = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
      dst->get( x, y, z, d.toIdx() ) = pdf;
      pdfs[ d.toIdx() ]              = pdf;
   }

   Vector3<real_t> velocity;
   real_t rho = dst->getDensityAndVelocity( velocity, x, y, z );

   // collide
   for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
   {
      const real_t fsym  = EquilibriumDistribution< LatticeModel_T >::getSymmetricPart ( *d, velocity, rho );
      const real_t fasym = EquilibriumDistribution< LatticeModel_T >::getAsymmetricPart( *d, velocity, rho );

      const real_t f     = pdfs[ d.toIdx() ];
      const real_t finv  = pdfs[ d.toInvIdx() ];

      dst->get( x, y, z, d.toIdx() ) = f - lambda_e_ * ( real_t( 0.5 ) * ( f + finv ) - fsym )
                                         - lambda_d_ * ( real_t( 0.5 ) * ( f - finv ) - fasym );
   }
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                    ! std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   real_t pdfs[ Stencil::Size ];
   
   // stream pull
   for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
   {
      const auto pdf = src.neighbor( d.inverseDir(), d.toIdx() );
      dst[ d.toIdx() ] = pdf;
      pdfs[ d.toIdx() ] = pdf;
   }

   Vector3<real_t> velocity;
   real_t rho = getDensityAndVelocity( velocity, *latticeModel_, dst );

   // collide
   for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
   {
      const real_t fsym  = EquilibriumDistribution< LatticeModel_T >::getSymmetricPart ( *d, velocity, rho );
      const real_t fasym = EquilibriumDistribution< LatticeModel_T >::getAsymmetricPart( *d, velocity, rho );

      const real_t f     = pdfs[ d.toIdx() ];
      const real_t finv  = pdfs[ d.toInvIdx() ];

      dst[ d.toIdx() ] = f - lambda_e_ * ( real_c( 0.5 ) * ( f + finv ) - fsym )
                           - lambda_d_ * ( real_c( 0.5 ) * ( f - finv ) - fasym );
   }
}



///////////////////////////////
// Specialization for D3Q19: //
// - incompressible          //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                     std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     ! LatticeModel_T::compressible &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value), "Only works with TRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() :
      lambda_e_( real_t(0) ), lambda_e_scaled_( real_t(0) ), lambda_d_scaled_( real_t(0) ),
      t0_( real_t(1.0) / real_t(3.0) ),
      t1x2_( real_t(1.0) / real_t(18.0) * real_t(2.0) ),
      t2x2_( real_t(1.0) / real_t(36.0) * real_t(2.0) ),
      fac1_( (real_t(1.0) / real_t(18.0) * real_t(2.0)) * (real_t(9.0) / real_t(2.0)) ),
      fac2_( (real_t(1.0) / real_t(36.0) * real_t(2.0)) * (real_t(9.0) / real_t(2.0)) ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      lambda_e_ = latticeModel.collisionModel().lambda_e();
      lambda_e_scaled_ = real_t(0.5) * lambda_e_;
      lambda_d_scaled_ = real_t(0.5) * latticeModel.collisionModel().lambda_d();
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;

private:

   real_t lambda_e_;
   real_t lambda_e_scaled_;
   real_t lambda_d_scaled_;

   const real_t t0_  ;
   const real_t t1x2_;
   const real_t t2x2_;
   const real_t fac1_;
   const real_t fac2_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
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

   const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

   dst->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e_) + lambda_e_ * t0_ * feq_common;

   const real_t velXPY = velX + velY;
   const real_t  sym_NE_SW = lambda_e_scaled_ * ( dd_tmp_NE + dd_tmp_SW - fac2_ * velXPY * velXPY - t2x2_ * feq_common );
   const real_t asym_NE_SW = lambda_d_scaled_ * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2_ * velXPY );
   dst->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
   dst->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;

   const real_t velXMY = velX - velY;
   const real_t  sym_SE_NW = lambda_e_scaled_ * ( dd_tmp_SE + dd_tmp_NW - fac2_ * velXMY * velXMY - t2x2_ * feq_common );
   const real_t asym_SE_NW = lambda_d_scaled_ * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2_ * velXMY );
   dst->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
   dst->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;

   const real_t velXPZ = velX + velZ;
   const real_t  sym_TE_BW = lambda_e_scaled_ * ( dd_tmp_TE + dd_tmp_BW - fac2_ * velXPZ * velXPZ - t2x2_ * feq_common );
   const real_t asym_TE_BW = lambda_d_scaled_ * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2_ * velXPZ );
   dst->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
   dst->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;

   const real_t velXMZ = velX - velZ;
   const real_t  sym_BE_TW = lambda_e_scaled_ * ( dd_tmp_BE + dd_tmp_TW - fac2_ * velXMZ * velXMZ - t2x2_ * feq_common );
   const real_t asym_BE_TW = lambda_d_scaled_ * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2_ * velXMZ );
   dst->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
   dst->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;

   const real_t velYPZ = velY + velZ;
   const real_t  sym_TN_BS = lambda_e_scaled_ * ( dd_tmp_TN + dd_tmp_BS - fac2_ * velYPZ * velYPZ - t2x2_ * feq_common );
   const real_t asym_TN_BS = lambda_d_scaled_ * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2_ * velYPZ );
   dst->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
   dst->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;

   const real_t velYMZ = velY - velZ;
   const real_t  sym_BN_TS = lambda_e_scaled_ * ( dd_tmp_BN + dd_tmp_TS - fac2_ * velYMZ * velYMZ - t2x2_ * feq_common );
   const real_t asym_BN_TS = lambda_d_scaled_ * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2_ * velYMZ );
   dst->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
   dst->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;

   const real_t  sym_N_S = lambda_e_scaled_ * ( dd_tmp_N + dd_tmp_S - fac1_ * velY * velY - t1x2_ * feq_common );
   const real_t asym_N_S = lambda_d_scaled_ * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2_ * velY );
   dst->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
   dst->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;

   const real_t  sym_E_W = lambda_e_scaled_ * ( dd_tmp_E + dd_tmp_W - fac1_ * velX * velX - t1x2_ * feq_common );
   const real_t asym_E_W = lambda_d_scaled_ * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2_ * velX );
   dst->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
   dst->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;

   const real_t  sym_T_B = lambda_e_scaled_ * ( dd_tmp_T + dd_tmp_B  - fac1_ * velZ * velZ - t1x2_ * feq_common );
   const real_t asym_T_B = lambda_d_scaled_ * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2_ * velZ );
   dst->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
   dst->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
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

   const real_t feq_common = rho - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

   dst[ Stencil::idx[C] ]= dd_tmp_C * (real_t(1.0) - lambda_e_) + lambda_e_ * t0_ * feq_common;

   const real_t velXPY = velX + velY;
   const real_t  sym_NE_SW = lambda_e_scaled_ * ( dd_tmp_NE + dd_tmp_SW - fac2_ * velXPY * velXPY - t2x2_ * feq_common );
   const real_t asym_NE_SW = lambda_d_scaled_ * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2_ * velXPY );
   dst[ Stencil::idx[NE] ]= dd_tmp_NE - sym_NE_SW - asym_NE_SW;
   dst[ Stencil::idx[SW] ]= dd_tmp_SW - sym_NE_SW + asym_NE_SW;

   const real_t velXMY = velX - velY;
   const real_t  sym_SE_NW = lambda_e_scaled_ * ( dd_tmp_SE + dd_tmp_NW - fac2_ * velXMY * velXMY - t2x2_ * feq_common );
   const real_t asym_SE_NW = lambda_d_scaled_ * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2_ * velXMY );
   dst[ Stencil::idx[SE] ]= dd_tmp_SE - sym_SE_NW - asym_SE_NW;
   dst[ Stencil::idx[NW] ]= dd_tmp_NW - sym_SE_NW + asym_SE_NW;

   const real_t velXPZ = velX + velZ;
   const real_t  sym_TE_BW = lambda_e_scaled_ * ( dd_tmp_TE + dd_tmp_BW - fac2_ * velXPZ * velXPZ - t2x2_ * feq_common );
   const real_t asym_TE_BW = lambda_d_scaled_ * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2_ * velXPZ );
   dst[ Stencil::idx[TE] ]= dd_tmp_TE - sym_TE_BW - asym_TE_BW;
   dst[ Stencil::idx[BW] ]= dd_tmp_BW - sym_TE_BW + asym_TE_BW;

   const real_t velXMZ = velX - velZ;
   const real_t  sym_BE_TW = lambda_e_scaled_ * ( dd_tmp_BE + dd_tmp_TW - fac2_ * velXMZ * velXMZ - t2x2_ * feq_common );
   const real_t asym_BE_TW = lambda_d_scaled_ * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2_ * velXMZ );
   dst[ Stencil::idx[BE] ]= dd_tmp_BE - sym_BE_TW - asym_BE_TW;
   dst[ Stencil::idx[TW] ]= dd_tmp_TW - sym_BE_TW + asym_BE_TW;

   const real_t velYPZ = velY + velZ;
   const real_t  sym_TN_BS = lambda_e_scaled_ * ( dd_tmp_TN + dd_tmp_BS - fac2_ * velYPZ * velYPZ - t2x2_ * feq_common );
   const real_t asym_TN_BS = lambda_d_scaled_ * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2_ * velYPZ );
   dst[ Stencil::idx[TN] ]= dd_tmp_TN - sym_TN_BS - asym_TN_BS;
   dst[ Stencil::idx[BS] ]= dd_tmp_BS - sym_TN_BS + asym_TN_BS;

   const real_t velYMZ = velY - velZ;
   const real_t  sym_BN_TS = lambda_e_scaled_ * ( dd_tmp_BN + dd_tmp_TS - fac2_ * velYMZ * velYMZ - t2x2_ * feq_common );
   const real_t asym_BN_TS = lambda_d_scaled_ * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2_ * velYMZ );
   dst[ Stencil::idx[BN] ]= dd_tmp_BN - sym_BN_TS - asym_BN_TS;
   dst[ Stencil::idx[TS] ]= dd_tmp_TS - sym_BN_TS + asym_BN_TS;

   const real_t  sym_N_S = lambda_e_scaled_ * ( dd_tmp_N + dd_tmp_S - fac1_ * velY * velY - t1x2_ * feq_common );
   const real_t asym_N_S = lambda_d_scaled_ * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2_ * velY );
   dst[ Stencil::idx[N] ]= dd_tmp_N - sym_N_S - asym_N_S;
   dst[ Stencil::idx[S] ]= dd_tmp_S - sym_N_S + asym_N_S;

   const real_t  sym_E_W = lambda_e_scaled_ * ( dd_tmp_E + dd_tmp_W - fac1_ * velX * velX - t1x2_ * feq_common );
   const real_t asym_E_W = lambda_d_scaled_ * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2_ * velX );
   dst[ Stencil::idx[E] ]= dd_tmp_E - sym_E_W - asym_E_W;
   dst[ Stencil::idx[W] ]= dd_tmp_W - sym_E_W + asym_E_W;

   const real_t  sym_T_B = lambda_e_scaled_ * ( dd_tmp_T + dd_tmp_B  - fac1_ * velZ * velZ - t1x2_ * feq_common );
   const real_t asym_T_B = lambda_d_scaled_ * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2_ * velZ );
   dst[ Stencil::idx[T] ]= dd_tmp_T - sym_T_B - asym_T_B;
   dst[ Stencil::idx[B] ]= dd_tmp_B - sym_T_B + asym_T_B;
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                     std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                     LatticeModel_T::compressible &&
                                                                     std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                     >::type >
{
public:

   static_assert( (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value), "Only works with TRT!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() :
      lambda_e_( real_t(0) ), lambda_e_scaled_( real_t(0) ), lambda_d_scaled_( real_t(0) ),
      t0_0_( real_t(1.0) / real_t(3.0) ),
      t1x2_0_( real_t(1.0) / real_t(18.0) * real_t(2.0) ),
      t2x2_0_( real_t(1.0) / real_t(36.0) * real_t(2.0) ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      lambda_e_ = latticeModel.collisionModel().lambda_e();
      lambda_e_scaled_ = real_t(0.5) * lambda_e_;
      lambda_d_scaled_ = real_t(0.5) * latticeModel.collisionModel().lambda_d();
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;

private:

   real_t lambda_e_;
   real_t lambda_e_scaled_;
   real_t lambda_d_scaled_;

   const real_t t0_0_  ;
   const real_t t1x2_0_;
   const real_t t2x2_0_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
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

   const real_t feq_common = real_t(1.0) - real_t(1.5) * ( velX * velX + velY * velY + velZ * velZ );

   dst->get( x, y, z, Stencil::idx[C] ) = dd_tmp_C * (real_t(1.0) - lambda_e_) + lambda_e_ * t0_0_ * rho * feq_common;

   const real_t t2x2 = t2x2_0_ * rho;
   const real_t fac2 = t2x2 * (real_t(9.0) / real_t(2.0));

   const real_t velXPY = velX + velY;
   const real_t  sym_NE_SW = lambda_e_scaled_ * ( dd_tmp_NE + dd_tmp_SW - fac2 * velXPY * velXPY - t2x2 * feq_common );
   const real_t asym_NE_SW = lambda_d_scaled_ * ( dd_tmp_NE - dd_tmp_SW - real_t(3.0) * t2x2 * velXPY );
   dst->get( x, y, z, Stencil::idx[NE] ) = dd_tmp_NE - sym_NE_SW - asym_NE_SW;
   dst->get( x, y, z, Stencil::idx[SW] ) = dd_tmp_SW - sym_NE_SW + asym_NE_SW;

   const real_t velXMY = velX - velY;
   const real_t  sym_SE_NW = lambda_e_scaled_ * ( dd_tmp_SE + dd_tmp_NW - fac2 * velXMY * velXMY - t2x2 * feq_common );
   const real_t asym_SE_NW = lambda_d_scaled_ * ( dd_tmp_SE - dd_tmp_NW - real_t(3.0) * t2x2 * velXMY );
   dst->get( x, y, z, Stencil::idx[SE] ) = dd_tmp_SE - sym_SE_NW - asym_SE_NW;
   dst->get( x, y, z, Stencil::idx[NW] ) = dd_tmp_NW - sym_SE_NW + asym_SE_NW;

   const real_t velXPZ = velX + velZ;
   const real_t  sym_TE_BW = lambda_e_scaled_ * ( dd_tmp_TE + dd_tmp_BW - fac2 * velXPZ * velXPZ - t2x2 * feq_common );
   const real_t asym_TE_BW = lambda_d_scaled_ * ( dd_tmp_TE - dd_tmp_BW - real_t(3.0) * t2x2 * velXPZ );
   dst->get( x, y, z, Stencil::idx[TE] ) = dd_tmp_TE - sym_TE_BW - asym_TE_BW;
   dst->get( x, y, z, Stencil::idx[BW] ) = dd_tmp_BW - sym_TE_BW + asym_TE_BW;

   const real_t velXMZ = velX - velZ;
   const real_t  sym_BE_TW = lambda_e_scaled_ * ( dd_tmp_BE + dd_tmp_TW - fac2 * velXMZ * velXMZ - t2x2 * feq_common );
   const real_t asym_BE_TW = lambda_d_scaled_ * ( dd_tmp_BE - dd_tmp_TW - real_t(3.0) * t2x2 * velXMZ );
   dst->get( x, y, z, Stencil::idx[BE] ) = dd_tmp_BE - sym_BE_TW - asym_BE_TW;
   dst->get( x, y, z, Stencil::idx[TW] ) = dd_tmp_TW - sym_BE_TW + asym_BE_TW;

   const real_t velYPZ = velY + velZ;
   const real_t  sym_TN_BS = lambda_e_scaled_ * ( dd_tmp_TN + dd_tmp_BS - fac2 * velYPZ * velYPZ - t2x2 * feq_common );
   const real_t asym_TN_BS = lambda_d_scaled_ * ( dd_tmp_TN - dd_tmp_BS - real_t(3.0) * t2x2 * velYPZ );
   dst->get( x, y, z, Stencil::idx[TN] ) = dd_tmp_TN - sym_TN_BS - asym_TN_BS;
   dst->get( x, y, z, Stencil::idx[BS] ) = dd_tmp_BS - sym_TN_BS + asym_TN_BS;

   const real_t velYMZ = velY - velZ;
   const real_t  sym_BN_TS = lambda_e_scaled_ * ( dd_tmp_BN + dd_tmp_TS - fac2 * velYMZ * velYMZ - t2x2 * feq_common );
   const real_t asym_BN_TS = lambda_d_scaled_ * ( dd_tmp_BN - dd_tmp_TS - real_t(3.0) * t2x2 * velYMZ );
   dst->get( x, y, z, Stencil::idx[BN] ) = dd_tmp_BN - sym_BN_TS - asym_BN_TS;
   dst->get( x, y, z, Stencil::idx[TS] ) = dd_tmp_TS - sym_BN_TS + asym_BN_TS;

   const real_t t1x2 = t1x2_0_ * rho;
   const real_t fac1 = t1x2 * (real_t(9.0) / real_t(2.0));

   const real_t  sym_N_S = lambda_e_scaled_ * ( dd_tmp_N + dd_tmp_S - fac1 * velY * velY - t1x2 * feq_common );
   const real_t asym_N_S = lambda_d_scaled_ * ( dd_tmp_N - dd_tmp_S - real_t(3.0) * t1x2 * velY );
   dst->get( x, y, z, Stencil::idx[N] ) = dd_tmp_N - sym_N_S - asym_N_S;
   dst->get( x, y, z, Stencil::idx[S] ) = dd_tmp_S - sym_N_S + asym_N_S;

   const real_t  sym_E_W = lambda_e_scaled_ * ( dd_tmp_E + dd_tmp_W - fac1 * velX * velX - t1x2 * feq_common );
   const real_t asym_E_W = lambda_d_scaled_ * ( dd_tmp_E - dd_tmp_W - real_t(3.0) * t1x2 * velX );
   dst->get( x, y, z, Stencil::idx[E] ) = dd_tmp_E - sym_E_W - asym_E_W;
   dst->get( x, y, z, Stencil::idx[W] ) = dd_tmp_W - sym_E_W + asym_E_W;

   const real_t  sym_T_B = lambda_e_scaled_ * ( dd_tmp_T + dd_tmp_B  - fac1 * velZ * velZ - t1x2 * feq_common );
   const real_t asym_T_B = lambda_d_scaled_ * ( dd_tmp_T - dd_tmp_B - real_t(3.0) * t1x2 * velZ );
   dst->get( x, y, z, Stencil::idx[T] ) = dd_tmp_T - sym_T_B - asym_T_B;
   dst->get( x, y, z, Stencil::idx[B] ) = dd_tmp_B - sym_T_B + asym_T_B;
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value &&
                                                                    std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value &&
                                                                    LatticeModel_T::compressible &&
                                                                    std::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value
                                                                    >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   operator()( const_cast< PdfField_T * >( static_cast< const PdfField_T * >(src.getField()) ),
               const_cast< PdfField_T * >( static_cast< const PdfField_T * >(dst.getField()) ),
               src.x(), src.y(), src.z() );
}



} // namespace lbm
} // namespace walberla
