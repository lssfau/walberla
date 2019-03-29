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
//! \file AdvectionDiffusionSweep.impl.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
// @see Chopard et.al.: The lattice Boltzmann advection-diffusion model revisited
//
//======================================================================================================================

#include "core/OpenMP.h"
#include "domain_decomposition/BlockDataID.h"
#include "field/Field.h"
#include "field/FieldClone.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/sweeps/FlagFieldSweepBase.h"
#include "lbm/sweeps/Streaming.h"

#include <type_traits>


namespace walberla {
namespace lbm {

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T, class Enable = void >
class AdvectionDiffusionSweep
{
	static_assert(never_true<LM_AdvDiff>::value, "Instantiating 'lbm::AdvectionDiffusionSweep' failed");
};



template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  LM_AdvDiff::CollisionModel::constant &&
                                  LM_AdvDiff::compressible &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                  ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                               >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;
     
   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ){}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ){}

   void operator() ( IBlock * block );
   
   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );

protected:
   const ConstBlockDataID velID_;
};

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void  AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  LM_AdvDiff::CollisionModel::constant &&
                                  LM_AdvDiff::compressible &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                  ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                               >::type
   > ::operator() ( IBlock * block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, dst, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // stream & collide

   const auto & lm = src->latticeModel( );
   dst->resetLatticeModel( lm ); // required so that the member function 'getDensity' can be called for dst!

   const real_t omega_eq = lm.collisionModel( ).omega( );
   const real_t omega_st = real_t( 1 ) - omega_eq;

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // stream
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
                  dst->get( x, y, z, d.toIdx( ) ) = src->get( x - d.cx( ), y - d.cy( ), z - d.cz( ), d.toIdx( ) );

               // macroscopic
               const auto velocity = vel->get( x, y, z );
               const real_t scalar = dst->getDensity( x, y, z );

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  dst->get( x, y, z, d.toIdx( ) ) = omega_st * dst->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
               }
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flag );

   dst->resetLatticeModel( src->latticeModel( ) ); // required so that the member function 'getDensity' can be called for dst!

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flag, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T   * src( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   const real_t omega_eq = src->latticeModel( ).collisionModel( ).omega( );
   const real_t omega_st = real_t( 1 ) - omega_eq;

   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // macroscopic
               const auto velocity = vel->get( x, y, z );
               const real_t scalar = src->getDensity( x, y, z );

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  src->get( x, y, z, d.toIdx( ) ) = omega_st * src->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
               }
            }
         }
      }
   }
}



template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  ! LM_AdvDiff::CollisionModel::constant &&
                                  LM_AdvDiff::compressible &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                  ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                               >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;
      
   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ){}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ){}
   
   void operator() ( IBlock * block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );

protected:
   const ConstBlockDataID velID_;
};


template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 !( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::operator() ( IBlock * block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, dst, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // stream & collide

   const auto & lm = src->latticeModel( );
   dst->resetLatticeModel( lm ); // required so that the member function 'getDensity' can be called for dst!

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // stream
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
                  dst->get( x, y, z, d.toIdx( ) ) = src->get( x - d.cx( ), y - d.cy( ), z - d.cz( ), d.toIdx( ) );

               // macroscopic
               const auto velocity = vel->get( x, y, z );
               const real_t scalar = dst->getDensity( x, y, z );

               const real_t omega_eq = lm.collisionModel( ).omega( x, y, z );
               const real_t omega_st = real_t( 1 ) - omega_eq;

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  dst->get( x, y, z, d.toIdx( ) ) = omega_st * dst->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
               }
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flag );

   dst->resetLatticeModel( src->latticeModel( ) ); // required so that the member function 'getDensity' can be called for dst!

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flag, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T   * src( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // macroscopic
               const auto velocity = vel->get( x, y, z );
               const real_t scalar = src->getDensity( x, y, z );

               const real_t omega_eq = src->latticeModel( ).collisionModel( ).omega( x, y, z );
               const real_t omega_st = real_t( 1 ) - omega_eq;

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  src->get( x, y, z, d.toIdx( ) ) = omega_st * src->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
               }
            }
         }
      }
   }
}



template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  LM_AdvDiff::CollisionModel::constant &&
                                  LM_AdvDiff::compressible &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                  ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                               >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value),  "Only works with correction force!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;
     
   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask, const BlockDataID & oldMomDensity)
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ), oldMomDensity_( oldMomDensity ) {}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask, const BlockDataID & oldMomDensity)
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ), oldMomDensity_( oldMomDensity ) {}

   void operator() ( IBlock * block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   
protected:
   const ConstBlockDataID velID_;
   const BlockDataID oldMomDensity_;
};

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1)
                              >::type
   > ::operator() ( IBlock * block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, dst, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   VectorField_T * oldMomDensityField = block->getData<VectorField_T>( oldMomDensity_ );

   WALBERLA_ASSERT_NOT_NULLPTR( oldMomDensityField );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->zSize( ), src->zSize( ) );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

   // stream & collide
   const auto & lm = src->latticeModel( );
   dst->resetLatticeModel( lm ); // required so that the member function 'getDensity' can be called for dst!

   const real_t omega_eq = lm.collisionModel( ).omega( );
   const real_t omega_st = real_t( 1 ) - omega_eq;
   const real_t omega_dt = real_t( 3 ) - real_c( 1.5 )*omega_eq;

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // stream
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
                  dst->get( x, y, z, d.toIdx( ) ) = src->get( x - d.cx( ), y - d.cy( ), z - d.cz( ), d.toIdx( ) );

               // macroscopic
               const auto   tmp = vel->get( x, y, z );
               const real_t scl = dst->getDensity( x, y, z );

               const auto scl_vel = scl * tmp;
               typename VectorField_T::value_type & scl_vel_old = oldMomDensityField->get( x, y, z );

               const auto dt_scl_vel = omega_dt * ( scl_vel - scl_vel_old );
               scl_vel_old = scl_vel;

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  dst->get( x, y, z, d.toIdx( ) ) = omega_st * dst->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, tmp, scl ) +
                     LM_AdvDiff::w[d.toIdx( )] * ( d.cx( )*dt_scl_vel[0] + d.cy( )*dt_scl_vel[1] + d.cz( )*dt_scl_vel[2] );
               }
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flag( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flag );

   dst->resetLatticeModel( src->latticeModel( ) ); // required so that the member function 'getDensity' can be called for dst!

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flag, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 LM_AdvDiff::compressible &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 ! ( std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value && LM_AdvDiff::equilibriumAccuracyOrder == 1 )
                              >::type
   > ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T   * src( NULL );
   const FlagField_T * flag( NULL );

   auto activeFlag = this->getLbmMaskAndFields( block, src, flag );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flag );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   VectorField_T * oldMomDensityField = block->getData<VectorField_T>( oldMomDensity_ );

   WALBERLA_ASSERT_NOT_NULLPTR( oldMomDensityField );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->zSize( ), src->zSize( ) );

   const real_t omega_eq = src->latticeModel( ).collisionModel( ).omega( );
   const real_t omega_st = real_t( 1 ) - omega_eq;
   const real_t omega_dt = real_t( 3 ) - real_c( 1.5 )*omega_eq;

   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            if( flag->isPartOfMaskSet( x, y, z, activeFlag ) )
            {
               // macroscopic
               const auto   tmp = vel->get( x, y, z );
               const real_t scl = src->getDensity( x, y, z );

               const auto scl_vel = scl * tmp;
               typename VectorField_T::value_type & scl_vel_old = oldMomDensityField->get( x, y, z );

               const auto dt_scl_vel = omega_dt * ( scl_vel - scl_vel_old );
               scl_vel_old = scl_vel;

               // collide
               for( auto d = Stencil::begin( ); d != Stencil::end( ); ++d )
               {
                  src->get( x, y, z, d.toIdx( ) ) = omega_st * src->get( x, y, z, d.toIdx( ) ) +
                     omega_eq * EquilibriumDistribution< LM_AdvDiff >::get( *d, tmp, scl ) +
                     LM_AdvDiff::w[d.toIdx( )] * ( d.cx( )*dt_scl_vel[0] + d.cy( )*dt_scl_vel[1] + d.cz( )*dt_scl_vel[2] );
               }
            }
         }
      }
   }
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  LM_AdvDiff::CollisionModel::constant &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                  std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19>::value &&
                                  LM_AdvDiff::compressible &&
                                  LM_AdvDiff::equilibriumAccuracyOrder == 1
                               >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LM_AdvDiff::equilibriumAccuracyOrder == 1,                                                     "Only works with equilibrium order 1!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;
     
   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ){}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ){}
   
   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );

protected:
   const ConstBlockDataID velID_;
};

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19>::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type
   > ::operator()( IBlock * const block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // constants used during stream/collide

   const real_t omega = src->latticeModel( ).collisionModel( ).omega( );

   const real_t omega_trm( real_t( 1 ) - omega );
   const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
   const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
   const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );
   const real_t one_third( real_t( 1 ) / real_t( 3 ) );

   // stream & collide

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x - 1, y - 1, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y - 1, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x + 1, y - 1, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x + 1, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x + 1, y + 1, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y + 1, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x - 1, y + 1, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x - 1, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z - 1, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x - 1, y, z - 1, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y - 1, z - 1, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x + 1, y, z - 1, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y + 1, z - 1, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z + 1, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x - 1, y, z + 1, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y - 1, z + 1, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x + 1, y, z + 1, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y + 1, z + 1, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               dst->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t omega_w1_rho_velX = omega_w1_rho * velX;
               const real_t omega_w1_rho_velY = omega_w1_rho * velY;
               const real_t omega_w1_rho_velZ = omega_w1_rho * velZ;

               dst->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + omega_w1_rho_velX;
               dst->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - omega_w1_rho_velX;
               dst->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + omega_w1_rho_velY;
               dst->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - omega_w1_rho_velY;
               dst->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + omega_w1_rho_velZ;
               dst->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - omega_w1_rho_velZ;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t omega_w2_rho_velXmY = omega_w2_rho * ( velX - velY );

               dst->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - omega_w2_rho_velXmY;
               dst->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + omega_w2_rho_velXmY;

               const real_t omega_w2_rho_velXpY = omega_w2_rho * ( velX + velY );

               dst->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + omega_w2_rho_velXpY;
               dst->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - omega_w2_rho_velXpY;

               const real_t omega_w2_rho_velXmZ = omega_w2_rho * ( velX - velZ );

               dst->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - omega_w2_rho_velXmZ;
               dst->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + omega_w2_rho_velXmZ;

               const real_t omega_w2_rho_velXpZ = omega_w2_rho * ( velX + velZ );

               dst->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + omega_w2_rho_velXpZ;
               dst->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - omega_w2_rho_velXpZ;

               const real_t omega_w2_rho_velYmZ = omega_w2_rho * ( velY - velZ );

               dst->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - omega_w2_rho_velYmZ;
               dst->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + omega_w2_rho_velYmZ;

               const real_t omega_w2_rho_velYpZ = omega_w2_rho * ( velY + velZ );

               dst->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + omega_w2_rho_velYpZ;
               dst->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - omega_w2_rho_velYpZ;
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19>::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type
> ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flagField, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19>::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type
> ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers( ), numberOfGhostLayersToInclude );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // constants used during collide

   const real_t omega = src->latticeModel( ).collisionModel( ).omega( );

   const real_t omega_trm( real_t( 1 ) - omega );
   const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
   const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
   const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );
   const real_t one_third( real_t( 1 ) / real_t( 3 ) );

   // collide
   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               src->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t omega_w1_rho_velX = omega_w1_rho * velX;
               const real_t omega_w1_rho_velY = omega_w1_rho * velY;
               const real_t omega_w1_rho_velZ = omega_w1_rho * velZ;

               src->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + omega_w1_rho_velX;
               src->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - omega_w1_rho_velX;
               src->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + omega_w1_rho_velY;
               src->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - omega_w1_rho_velY;
               src->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + omega_w1_rho_velZ;
               src->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - omega_w1_rho_velZ;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t omega_w2_rho_velXmY = omega_w2_rho * ( velX - velY );

               src->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - omega_w2_rho_velXmY;
               src->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + omega_w2_rho_velXmY;

               const real_t omega_w2_rho_velXpY = omega_w2_rho * ( velX + velY );

               src->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + omega_w2_rho_velXpY;
               src->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - omega_w2_rho_velXpY;

               const real_t omega_w2_rho_velXmZ = omega_w2_rho * ( velX - velZ );

               src->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - omega_w2_rho_velXmZ;
               src->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + omega_w2_rho_velXmZ;

               const real_t omega_w2_rho_velXpZ = omega_w2_rho * ( velX + velZ );

               src->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + omega_w2_rho_velXpZ;
               src->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - omega_w2_rho_velXpZ;

               const real_t omega_w2_rho_velYmZ = omega_w2_rho * ( velY - velZ );

               src->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - omega_w2_rho_velYmZ;
               src->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + omega_w2_rho_velYmZ;

               const real_t omega_w2_rho_velYpZ = omega_w2_rho * ( velY + velZ );

               src->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + omega_w2_rho_velYpZ;
               src->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - omega_w2_rho_velYpZ;
            }
         }
      }
   }
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                               typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                  ! LM_AdvDiff::CollisionModel::constant &&
                                  std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                  std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                  LM_AdvDiff::compressible &&
                                  LM_AdvDiff::equilibriumAccuracyOrder == 1
                               >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LM_AdvDiff::equilibriumAccuracyOrder == 1,                                                     "Only works with equilibrium order 1!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;
      
   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ){}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask )
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ){}
   
   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );

protected:
   const ConstBlockDataID velID_;
};

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::operator()( IBlock * const block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // stream & collide

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               // constants used during stream/collide
               const real_t omega = src->latticeModel( ).collisionModel( ).omega( x, y, z );

               const real_t omega_trm( real_t( 1 ) - omega );
               const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
               const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
               const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );
               const real_t one_third( real_t( 1 ) / real_t( 3 ) );

               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x - 1, y - 1, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y - 1, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x + 1, y - 1, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x + 1, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x + 1, y + 1, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y + 1, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x - 1, y + 1, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x - 1, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z - 1, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x - 1, y, z - 1, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y - 1, z - 1, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x + 1, y, z - 1, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y + 1, z - 1, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z + 1, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x - 1, y, z + 1, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y - 1, z + 1, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x + 1, y, z + 1, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y + 1, z + 1, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               dst->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t omega_w1_rho_velX = omega_w1_rho * velX;
               const real_t omega_w1_rho_velY = omega_w1_rho * velY;
               const real_t omega_w1_rho_velZ = omega_w1_rho * velZ;

               dst->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + omega_w1_rho_velX;
               dst->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - omega_w1_rho_velX;
               dst->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + omega_w1_rho_velY;
               dst->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - omega_w1_rho_velY;
               dst->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + omega_w1_rho_velZ;
               dst->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - omega_w1_rho_velZ;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t omega_w2_rho_velXmY = omega_w2_rho * ( velX - velY );

               dst->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - omega_w2_rho_velXmY;
               dst->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + omega_w2_rho_velXmY;

               const real_t omega_w2_rho_velXpY = omega_w2_rho * ( velX + velY );

               dst->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + omega_w2_rho_velXpY;
               dst->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - omega_w2_rho_velXpY;

               const real_t omega_w2_rho_velXmZ = omega_w2_rho * ( velX - velZ );

               dst->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - omega_w2_rho_velXmZ;
               dst->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + omega_w2_rho_velXmZ;

               const real_t omega_w2_rho_velXpZ = omega_w2_rho * ( velX + velZ );

               dst->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + omega_w2_rho_velXpZ;
               dst->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - omega_w2_rho_velXpZ;

               const real_t omega_w2_rho_velYmZ = omega_w2_rho * ( velY - velZ );

               dst->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - omega_w2_rho_velYmZ;
               dst->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + omega_w2_rho_velYmZ;

               const real_t omega_w2_rho_velYpZ = omega_w2_rho * ( velY + velZ );

               dst->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + omega_w2_rho_velYpZ;
               dst->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - omega_w2_rho_velYpZ;
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flagField, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 ! LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers( ), numberOfGhostLayersToInclude );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   // collide
   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               // constants used during collide
               const real_t omega = src->latticeModel( ).collisionModel( ).omega( x, y, z );

               const real_t omega_trm( real_t( 1 ) - omega );
               const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
               const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
               const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );
               const real_t one_third( real_t( 1 ) / real_t( 3 ) );

               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               src->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t omega_w1_rho_velX = omega_w1_rho * velX;
               const real_t omega_w1_rho_velY = omega_w1_rho * velY;
               const real_t omega_w1_rho_velZ = omega_w1_rho * velZ;

               src->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + omega_w1_rho_velX;
               src->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - omega_w1_rho_velX;
               src->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + omega_w1_rho_velY;
               src->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - omega_w1_rho_velY;
               src->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + omega_w1_rho_velZ;
               src->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - omega_w1_rho_velZ;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t omega_w2_rho_velXmY = omega_w2_rho * ( velX - velY );

               src->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - omega_w2_rho_velXmY;
               src->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + omega_w2_rho_velXmY;

               const real_t omega_w2_rho_velXpY = omega_w2_rho * ( velX + velY );

               src->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + omega_w2_rho_velXpY;
               src->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - omega_w2_rho_velXpY;

               const real_t omega_w2_rho_velXmZ = omega_w2_rho * ( velX - velZ );

               src->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - omega_w2_rho_velXmZ;
               src->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + omega_w2_rho_velXmZ;

               const real_t omega_w2_rho_velXpZ = omega_w2_rho * ( velX + velZ );

               src->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + omega_w2_rho_velXpZ;
               src->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - omega_w2_rho_velXpZ;

               const real_t omega_w2_rho_velYmZ = omega_w2_rho * ( velY - velZ );

               src->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - omega_w2_rho_velYmZ;
               src->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + omega_w2_rho_velYmZ;

               const real_t omega_w2_rho_velYpZ = omega_w2_rho * ( velY + velZ );

               src->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + omega_w2_rho_velYpZ;
               src->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - omega_w2_rho_velYpZ;
            }
         }
      }
   }
}


///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - correction term         //
///////////////////////////////
template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
class AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
   typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                            LM_AdvDiff::CollisionModel::constant &&
                            std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                            std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                            LM_AdvDiff::compressible &&
                            LM_AdvDiff::equilibriumAccuracyOrder == 1
                            >::type > :
   public FlagFieldSweepBase< LM_AdvDiff, FlagField_T >
{
public:
   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value),  "Only works with correction force!" );
   static_assert( LM_AdvDiff::equilibriumAccuracyOrder == 1,                                                     "Only works with equilibrium order 1!" );

   typedef typename FlagFieldSweepBase<LM_AdvDiff,FlagField_T>::PdfField_T  AdvDiffPdfField_T;
   typedef typename LM_AdvDiff::Stencil                                     Stencil;


   AdvectionDiffusionSweep( const BlockDataID & advDiffID, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask, const BlockDataID & oldMomDensity)
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffID, flagID, lbmMask ), velID_( velID ), oldMomDensity_( oldMomDensity ) {}
   
   AdvectionDiffusionSweep( const BlockDataID & advDiffSrcID, const BlockDataID & advDiffDstIS, const ConstBlockDataID & velID, const ConstBlockDataID & flagID, const Set< FlagUID > & lbmMask, const BlockDataID & oldMomDensity)
      :  FlagFieldSweepBase<LM_AdvDiff,FlagField_T>( advDiffSrcID, advDiffDstIS, flagID, lbmMask ), velID_( velID ), oldMomDensity_( oldMomDensity ) {}

   void operator()( IBlock * const block );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t( 0u ) );

protected:
   const ConstBlockDataID velID_;
   const BlockDataID oldMomDensity_;
};

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::operator()( IBlock * const block )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), 1 );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   VectorField_T * oldMomDensityField = block->getData<VectorField_T>( oldMomDensity_ );

   WALBERLA_ASSERT_NOT_NULLPTR( oldMomDensityField );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->zSize( ), src->zSize( ) );

   // constants used during stream/collide

   const real_t omega = src->latticeModel( ).collisionModel( ).omega( );

   const real_t omega_trm( real_t( 1 ) - omega );
   const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
   const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
   const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );

   const real_t omega_dt1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * ( real_t( 1 ) - ( real_t( 1 ) / real_t( 2 ) ) * omega ) );
   const real_t omega_dt2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * ( real_t( 1 ) - ( real_t( 1 ) / real_t( 2 ) ) * omega ) );

   const real_t one_third( real_t( 1 ) / real_t( 3 ) );

   // stream & collide

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = 0; y < ySize; ++y ) {
         for( cell_idx_t x = 0; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x - 1, y - 1, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y - 1, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x + 1, y - 1, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x + 1, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x + 1, y + 1, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y + 1, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x - 1, y + 1, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x - 1, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z - 1, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x - 1, y, z - 1, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y - 1, z - 1, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x + 1, y, z - 1, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y + 1, z - 1, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z + 1, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x - 1, y, z + 1, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y - 1, z + 1, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x + 1, y, z + 1, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y + 1, z + 1, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               const real_t rho_velX = rho * velX;
               const real_t rho_velY = rho * velY;
               const real_t rho_velZ = rho * velZ;

               typename VectorField_T::value_type & old_rho_vel = oldMomDensityField->get( x, y, z );

               const real_t dt_rho_velX = rho_velX - old_rho_vel[0];
               const real_t dt_rho_velY = rho_velY - old_rho_vel[1];
               const real_t dt_rho_velZ = rho_velZ - old_rho_vel[2];

               old_rho_vel[0] = rho_velX;
               old_rho_vel[1] = rho_velY;
               old_rho_vel[2] = rho_velZ;

               dst->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t vel_trm_E_W = omega_w1_rho * velX + omega_dt1 * dt_rho_velX;
               const real_t vel_trm_N_S = omega_w1_rho * velY + omega_dt1 * dt_rho_velY;
               const real_t vel_trm_T_B = omega_w1_rho * velZ + omega_dt1 * dt_rho_velZ;

               dst->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + vel_trm_E_W;
               dst->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - vel_trm_E_W;
               dst->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + vel_trm_N_S;
               dst->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - vel_trm_N_S;
               dst->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + vel_trm_T_B;
               dst->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - vel_trm_T_B;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t vel_trm_NW_SE = omega_w2_rho * ( velX - velY ) + omega_dt2 * ( dt_rho_velX - dt_rho_velY );

               dst->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - vel_trm_NW_SE;
               dst->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + vel_trm_NW_SE;

               const real_t vel_trm_NE_SW = omega_w2_rho * ( velX + velY ) + omega_dt2 * ( dt_rho_velX + dt_rho_velY );

               dst->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + vel_trm_NE_SW;
               dst->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - vel_trm_NE_SW;

               const real_t vel_trm_TW_BE = omega_w2_rho * ( velX - velZ ) + omega_dt2 * ( dt_rho_velX - dt_rho_velZ );

               dst->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - vel_trm_TW_BE;
               dst->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + vel_trm_TW_BE;

               const real_t vel_trm_TE_BW = omega_w2_rho * ( velX + velZ ) + omega_dt2 * ( dt_rho_velX + dt_rho_velZ );

               dst->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + vel_trm_TE_BW;
               dst->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - vel_trm_TE_BW;

               const real_t vel_trm_TS_BN = omega_w2_rho * ( velY - velZ ) + omega_dt2 * ( dt_rho_velY - dt_rho_velZ );

               dst->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - vel_trm_TS_BN;
               dst->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + vel_trm_TS_BN;

               const real_t vel_trm_TN_BS = omega_w2_rho * ( velY + velZ ) + omega_dt2 * ( dt_rho_velY + dt_rho_velZ );

               dst->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + vel_trm_TN_BS;
               dst->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - vel_trm_TN_BS;
            }
         }
      }
   }

   src->swapDataPointers( dst );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   AdvDiffPdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   Stream< LM_AdvDiff, FlagField_T >::execute( src, dst, flagField, lbm, numberOfGhostLayersToInclude );
}

template< typename LM_AdvDiff, typename VelocityAdapter_T, typename FlagField_T, typename VectorField_T >
void AdvectionDiffusionSweep< LM_AdvDiff, VelocityAdapter_T, FlagField_T, VectorField_T,
                              typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                 LM_AdvDiff::CollisionModel::constant &&
                                 std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value &&
                                 std::is_same< typename LM_AdvDiff::Stencil, stencil::D3Q19 >::value &&
                                 LM_AdvDiff::compressible &&
                                 LM_AdvDiff::equilibriumAccuracyOrder == 1
                              >::type > ::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   AdvDiffPdfField_T * src( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers( ), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers( ), numberOfGhostLayersToInclude );

   const VelocityAdapter_T * vel = block->getData<VelocityAdapter_T>( velID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( vel );
   WALBERLA_ASSERT_EQUAL( vel->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( vel->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( vel->zSize( ), src->zSize( ) );

   VectorField_T * oldMomDensityField = block->getData<VectorField_T>( oldMomDensity_ );

   WALBERLA_ASSERT_NOT_NULLPTR( oldMomDensityField );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->xSize( ), src->xSize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->ySize( ), src->ySize( ) );
   WALBERLA_ASSERT_EQUAL( oldMomDensityField->zSize( ), src->zSize( ) );

   // constants used during collide

   const real_t omega = src->latticeModel( ).collisionModel( ).omega( );

   const real_t omega_trm( real_t( 1 ) - omega );
   const real_t  omega_w0( real_t( 3 ) * ( real_t( 1 ) / real_t( 3 ) ) * omega );
   const real_t  omega_w1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * omega );
   const real_t  omega_w2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * omega );

   const real_t omega_dt1( real_t( 3 ) * ( real_t( 1 ) / real_t( 18 ) ) * ( real_t( 1 ) - ( real_t( 1 ) / real_t( 2 ) ) * omega ) );
   const real_t omega_dt2( real_t( 3 ) * ( real_t( 1 ) / real_t( 36 ) ) * ( real_t( 1 ) - ( real_t( 1 ) / real_t( 2 ) ) * omega ) );

   const real_t one_third( real_t( 1 ) / real_t( 3 ) );

   // collide
   const cell_idx_t start = -cell_idx_c( numberOfGhostLayersToInclude );

   const cell_idx_t xSize = cell_idx_c( src->xSize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t ySize = cell_idx_c( src->ySize( ) + numberOfGhostLayersToInclude );
   const cell_idx_t zSize = cell_idx_c( src->zSize( ) + numberOfGhostLayersToInclude );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
#pragma omp parallel for schedule(static)
   for( int iz = int_c( start ); iz < izSize; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
#else
   for( cell_idx_t z = start; z < zSize; ++z ) {
#endif
      for( cell_idx_t y = start; y < ySize; ++y ) {
         for( cell_idx_t x = start; x < xSize; ++x )
         {
            using namespace stencil;

            if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
            {
               const auto   tmp = vel->get( x, y, z );
               const real_t velX = tmp[0];
               const real_t velY = tmp[1];
               const real_t velZ = tmp[2];

               const real_t dd_tmp_NE = src->get( x, y, z, Stencil::idx[NE] );
               const real_t dd_tmp_N = src->get( x, y, z, Stencil::idx[N] );
               const real_t dd_tmp_NW = src->get( x, y, z, Stencil::idx[NW] );
               const real_t dd_tmp_W = src->get( x, y, z, Stencil::idx[W] );
               const real_t dd_tmp_SW = src->get( x, y, z, Stencil::idx[SW] );
               const real_t dd_tmp_S = src->get( x, y, z, Stencil::idx[S] );
               const real_t dd_tmp_SE = src->get( x, y, z, Stencil::idx[SE] );
               const real_t dd_tmp_E = src->get( x, y, z, Stencil::idx[E] );
               const real_t dd_tmp_T = src->get( x, y, z, Stencil::idx[T] );
               const real_t dd_tmp_TE = src->get( x, y, z, Stencil::idx[TE] );
               const real_t dd_tmp_TN = src->get( x, y, z, Stencil::idx[TN] );
               const real_t dd_tmp_TW = src->get( x, y, z, Stencil::idx[TW] );
               const real_t dd_tmp_TS = src->get( x, y, z, Stencil::idx[TS] );
               const real_t dd_tmp_B = src->get( x, y, z, Stencil::idx[B] );
               const real_t dd_tmp_BE = src->get( x, y, z, Stencil::idx[BE] );
               const real_t dd_tmp_BN = src->get( x, y, z, Stencil::idx[BN] );
               const real_t dd_tmp_BW = src->get( x, y, z, Stencil::idx[BW] );
               const real_t dd_tmp_BS = src->get( x, y, z, Stencil::idx[BS] );
               const real_t dd_tmp_C = src->get( x, y, z, Stencil::idx[C] );

               const real_t rho = dd_tmp_C + dd_tmp_N + dd_tmp_S + dd_tmp_W + dd_tmp_E + dd_tmp_T + dd_tmp_B +
                  dd_tmp_NE + dd_tmp_NW + dd_tmp_SW + dd_tmp_SE + dd_tmp_TE + dd_tmp_TN + dd_tmp_TW + dd_tmp_TS + dd_tmp_BE + dd_tmp_BN + dd_tmp_BW + dd_tmp_BS;

               const real_t rho_velX = rho * velX;
               const real_t rho_velY = rho * velY;
               const real_t rho_velZ = rho * velZ;

               typename VectorField_T::value_type & old_rho_vel = oldMomDensityField->get( x, y, z );

               const real_t dt_rho_velX = rho_velX - old_rho_vel[0];
               const real_t dt_rho_velY = rho_velY - old_rho_vel[1];
               const real_t dt_rho_velZ = rho_velZ - old_rho_vel[2];

               old_rho_vel[0] = rho_velX;
               old_rho_vel[1] = rho_velY;
               old_rho_vel[2] = rho_velZ;

               src->get( x, y, z, Stencil::idx[C] ) = omega_trm * dd_tmp_C + omega_w0 * rho * one_third;

               const real_t omega_w1_rho = omega_w1     * rho;
               const real_t omega_w1_rho_third = omega_w1_rho * one_third;

               const real_t vel_trm_E_W = omega_w1_rho * velX + omega_dt1 * dt_rho_velX;
               const real_t vel_trm_N_S = omega_w1_rho * velY + omega_dt1 * dt_rho_velY;
               const real_t vel_trm_T_B = omega_w1_rho * velZ + omega_dt1 * dt_rho_velZ;

               src->get( x, y, z, Stencil::idx[E] ) = omega_trm * dd_tmp_E + omega_w1_rho_third + vel_trm_E_W;
               src->get( x, y, z, Stencil::idx[W] ) = omega_trm * dd_tmp_W + omega_w1_rho_third - vel_trm_E_W;
               src->get( x, y, z, Stencil::idx[N] ) = omega_trm * dd_tmp_N + omega_w1_rho_third + vel_trm_N_S;
               src->get( x, y, z, Stencil::idx[S] ) = omega_trm * dd_tmp_S + omega_w1_rho_third - vel_trm_N_S;
               src->get( x, y, z, Stencil::idx[T] ) = omega_trm * dd_tmp_T + omega_w1_rho_third + vel_trm_T_B;
               src->get( x, y, z, Stencil::idx[B] ) = omega_trm * dd_tmp_B + omega_w1_rho_third - vel_trm_T_B;

               const real_t omega_w2_rho = omega_w2     * rho;
               const real_t omega_w2_rho_third = omega_w2_rho * one_third;

               const real_t vel_trm_NW_SE = omega_w2_rho * ( velX - velY ) + omega_dt2 * ( dt_rho_velX - dt_rho_velY );

               src->get( x, y, z, Stencil::idx[NW] ) = omega_trm * dd_tmp_NW + omega_w2_rho_third - vel_trm_NW_SE;
               src->get( x, y, z, Stencil::idx[SE] ) = omega_trm * dd_tmp_SE + omega_w2_rho_third + vel_trm_NW_SE;

               const real_t vel_trm_NE_SW = omega_w2_rho * ( velX + velY ) + omega_dt2 * ( dt_rho_velX + dt_rho_velY );

               src->get( x, y, z, Stencil::idx[NE] ) = omega_trm * dd_tmp_NE + omega_w2_rho_third + vel_trm_NE_SW;
               src->get( x, y, z, Stencil::idx[SW] ) = omega_trm * dd_tmp_SW + omega_w2_rho_third - vel_trm_NE_SW;

               const real_t vel_trm_TW_BE = omega_w2_rho * ( velX - velZ ) + omega_dt2 * ( dt_rho_velX - dt_rho_velZ );

               src->get( x, y, z, Stencil::idx[TW] ) = omega_trm * dd_tmp_TW + omega_w2_rho_third - vel_trm_TW_BE;
               src->get( x, y, z, Stencil::idx[BE] ) = omega_trm * dd_tmp_BE + omega_w2_rho_third + vel_trm_TW_BE;

               const real_t vel_trm_TE_BW = omega_w2_rho * ( velX + velZ ) + omega_dt2 * ( dt_rho_velX + dt_rho_velZ );

               src->get( x, y, z, Stencil::idx[TE] ) = omega_trm * dd_tmp_TE + omega_w2_rho_third + vel_trm_TE_BW;
               src->get( x, y, z, Stencil::idx[BW] ) = omega_trm * dd_tmp_BW + omega_w2_rho_third - vel_trm_TE_BW;

               const real_t vel_trm_TS_BN = omega_w2_rho * ( velY - velZ ) + omega_dt2 * ( dt_rho_velY - dt_rho_velZ );

               src->get( x, y, z, Stencil::idx[TS] ) = omega_trm * dd_tmp_TS + omega_w2_rho_third - vel_trm_TS_BN;
               src->get( x, y, z, Stencil::idx[BN] ) = omega_trm * dd_tmp_BN + omega_w2_rho_third + vel_trm_TS_BN;

               const real_t vel_trm_TN_BS = omega_w2_rho * ( velY + velZ ) + omega_dt2 * ( dt_rho_velY + dt_rho_velZ );

               src->get( x, y, z, Stencil::idx[TN] ) = omega_trm * dd_tmp_TN + omega_w2_rho_third + vel_trm_TN_BS;
               src->get( x, y, z, Stencil::idx[BS] ) = omega_trm * dd_tmp_BS + omega_w2_rho_third - vel_trm_TN_BS;
            }
         }
      }
   }
}


} // namespace lbm
} // namespace walberla
