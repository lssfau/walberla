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
//! \file CGFixedStencilIteration.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Set.h"
#include "core/logging/Logging.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/BlockStorage.h"

#include "field/GhostLayerField.h"
#include "field/iterators/IteratorMacros.h"

#include <functional>



namespace walberla {
namespace pde {



template< typename Stencil_T >
class CGFixedStencilIteration
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   CGFixedStencilIteration( BlockStorage & blocks,
                            const BlockDataID & uId, const BlockDataID & rId, const BlockDataID & dId, const BlockDataID & zId,
                            const BlockDataID & fId, const std::vector< real_t > & weights,
                            const uint_t iterations, const std::function< void () > & synchronizeD,
                            const real_t residualNormThreshold = real_t(0),
                            const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
      
   void operator()();
   
protected:

   ////////////////////////////
   // building blocks for CG //
   ////////////////////////////
   void   calcR();                       // r = f - Au
   real_t scalarProductRR();             // r*r
   void   copyRToD();                    // d = r
   void   calcAd();                      // z = Ad
   real_t scalarProductDZ();             // d*z
   void   updateU( const real_t alpha ); // u = u + alpha * d
   void   updateR( const real_t alpha ); // r = r - alpha * z
   void   updateD( const real_t beta  ); // d = r + beta * d



   BlockStorage & blocks_;

   const BlockDataID uId_;
   const BlockDataID rId_;
   const BlockDataID dId_;
   const BlockDataID zId_;
   const BlockDataID fId_;
   
   real_t cells_;

   real_t w_[ Stencil_T::Size ];

   uint_t iterations_;
   real_t residualNormThreshold_;
   
   std::function< void () > synchronizeD_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



template< typename Stencil_T >
CGFixedStencilIteration< Stencil_T >::CGFixedStencilIteration( BlockStorage & blocks,
                                                              const BlockDataID & uId, const BlockDataID & rId, const BlockDataID & dId, const BlockDataID & zId,
                                                              const BlockDataID & fId, const std::vector< real_t > & weights,
                                                              const uint_t iterations, const std::function< void () > & synchronizeD,
                                                              const real_t residualNormThreshold,
                                                              const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors ) :
   blocks_( blocks ), uId_( uId ), rId_( rId ), dId_( dId ), zId_( zId ), fId_( fId ),
   iterations_( iterations ),
   residualNormThreshold_( residualNormThreshold ),
   synchronizeD_( synchronizeD ),
   requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
{
   WALBERLA_ASSERT_EQUAL( weights.size(), Stencil_T::Size );
   for( uint_t i = uint_t(0); i < Stencil_T::Size; ++i )
      w_[i] = weights[i];
      
   uint_t cells( uint_t(0) );
   
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      const Field_T * const u = block->template getData< const Field_T >( uId_ );
      cells += u->xyzSize().numCells();
   }
   
   cells_ = real_c( cells );
   mpi::allReduceInplace( cells_, mpi::SUM );         
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::operator()()
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting CG iteration with a maximum number of " << iterations_ << " iterations" );

   calcR(); // r = f - Au

   real_t rr0 = scalarProductRR(); // r*r
   real_t residualNorm = std::sqrt( rr0 / cells_ );
   
   if( residualNorm >= residualNormThreshold_ )
   {
      copyRToD(); // d = r
      
      uint_t i( uint_t(0) );
      while( i < iterations_ )
      {
         synchronizeD_();
         
         calcAd(); // z = Ad
         
         const real_t alpha = rr0 / scalarProductDZ(); // alpha = r*r / d*z
         updateU( alpha ); // u = u + alpha * d
         updateR( alpha ); // r = r - alpha * z
         
         const real_t rr1 = scalarProductRR();
         residualNorm = std::sqrt( rr1 / cells_ );
         if( residualNorm < residualNormThreshold_ )
         {
            WALBERLA_LOG_PROGRESS_ON_ROOT( "Aborting CG iteration (residual norm threshold reached):"
                                           "\n  residual norm threshold: " << residualNormThreshold_ <<
                                           "\n  residual norm:           " << residualNorm );
            break;
         }
         
         const real_t beta = rr1 / rr0; // beta = r*r (current) / r*r (previous)
         updateD( beta ); // d = r + beta * d
         
         rr0 = rr1;
         
         ++i;
      }
      
      WALBERLA_LOG_PROGRESS_ON_ROOT( "CG iteration finished after " << i << " iterations" );
   }
   else
   {
      WALBERLA_LOG_PROGRESS_ON_ROOT( "Aborting CG without a single iteration (residual norm threshold already reached):"
                                     "\n  residual norm threshold: " << residualNormThreshold_ <<
                                     "\n  residual norm:           " << residualNorm );   
   }   
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::calcR()  // r = f - Au
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * rf = block->template getData< Field_T >( rId_ );
      Field_T * ff = block->template getData< Field_T >( fId_ );
      Field_T * uf = block->template getData< Field_T >( uId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( rf );
      WALBERLA_ASSERT_NOT_NULLPTR( ff );
      WALBERLA_ASSERT_NOT_NULLPTR( uf );

      WALBERLA_ASSERT_EQUAL( rf->xyzSize(), ff->xyzSize() );
      WALBERLA_ASSERT_EQUAL( rf->xyzSize(), uf->xyzSize() );

      WALBERLA_ASSERT_GREATER_EQUAL( uf->nrOfGhostLayers(), 1 );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( uf, 

         rf->get(x,y,z) = ff->get(x,y,z);

         for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
            rf->get(x,y,z) -= w_[ dir.toIdx() ] * uf->getNeighbor( x, y, z, *dir );
      )
   }
}



template< typename Stencil_T >
real_t CGFixedStencilIteration< Stencil_T >::scalarProductRR() // r*r
{
   real_t result( real_t(0) );
   
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * rf = block->template getData< Field_T >( rId_ );
      
      real_t blockResult( real_t(0) );
      
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP( rf, omp parallel for schedule(static) reduction(+:blockResult),

         const real_t v = rf->get(x,y,z);
         blockResult += v * v;
      )
      
      result += blockResult;
   }
   
   mpi::allReduceInplace( result, mpi::SUM );
   return result;
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::copyRToD() // d = r
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * rf = block->template getData< Field_T >( rId_ );
      Field_T * df = block->template getData< Field_T >( dId_ );
      
      WALBERLA_ASSERT_NOT_NULLPTR( rf );
      WALBERLA_ASSERT_NOT_NULLPTR( df );

      WALBERLA_ASSERT_EQUAL( rf->xyzSize(), df->xyzSize() );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( rf, 

         df->get(x,y,z) = rf->get(x,y,z);
      )
   }
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::calcAd() // z = Ad
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * zf = block->template getData< Field_T >( zId_ );
      Field_T * df = block->template getData< Field_T >( dId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( zf );
      WALBERLA_ASSERT_NOT_NULLPTR( df );

      WALBERLA_ASSERT_EQUAL( zf->xyzSize(), df->xyzSize() );
      
      WALBERLA_ASSERT_GREATER_EQUAL( df->nrOfGhostLayers(), 1 );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( df,

         zf->get(x,y,z) = w_[ Stencil_T::idx[stencil::C] ] * df->get(x,y,z);
         
         for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
            zf->get(x,y,z) += w_[ dir.toIdx() ] * df->getNeighbor( x, y, z, *dir );
      )
   }
}



template< typename Stencil_T >
real_t CGFixedStencilIteration< Stencil_T >::scalarProductDZ() // d*z
{
   real_t result( real_t(0) );
   
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * df = block->template getData< Field_T >( dId_ );
      Field_T * zf = block->template getData< Field_T >( zId_ );
      
      WALBERLA_ASSERT_NOT_NULLPTR( df );
      WALBERLA_ASSERT_NOT_NULLPTR( zf );

      WALBERLA_ASSERT_EQUAL( df->xyzSize(), zf->xyzSize() );
      
      real_t blockResult( real_t(0) );
      
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP( df, omp parallel for schedule(static) reduction(+:blockResult),

         blockResult += df->get(x,y,z) * zf->get(x,y,z);
      )
      
      result += blockResult;
   }
   
   mpi::allReduceInplace( result, mpi::SUM );
   return result;
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::updateU( const real_t alpha ) // u = u + alpha * d
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * uf = block->template getData< Field_T >( uId_ );
      Field_T * df = block->template getData< Field_T >( dId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( uf );
      WALBERLA_ASSERT_NOT_NULLPTR( df );

      WALBERLA_ASSERT_EQUAL( uf->xyzSize(), df->xyzSize() );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( uf,

         uf->get(x,y,z) = uf->get(x,y,z) + alpha * df->get(x,y,z);
      )
   }
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::updateR( const real_t alpha ) // r = r - alpha * z
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * rf = block->template getData< Field_T >( rId_ );
      Field_T * zf = block->template getData< Field_T >( zId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( rf );
      WALBERLA_ASSERT_NOT_NULLPTR( zf );

      WALBERLA_ASSERT_EQUAL( rf->xyzSize(), zf->xyzSize() );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( rf,

         rf->get(x,y,z) = rf->get(x,y,z) - alpha * zf->get(x,y,z);
      )
   }
}



template< typename Stencil_T >
void CGFixedStencilIteration< Stencil_T >::updateD( const real_t beta  ) // d = r + beta * d
{
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      Field_T * df = block->template getData< Field_T >( dId_ );
      Field_T * rf = block->template getData< Field_T >( rId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( df );
      WALBERLA_ASSERT_NOT_NULLPTR( rf );

      WALBERLA_ASSERT_EQUAL( df->xyzSize(), rf->xyzSize() );
      
      WALBERLA_FOR_ALL_CELLS_XYZ( df,

         df->get(x,y,z) = rf->get(x,y,z) + beta * df->get(x,y,z);
      )
   }
}



} // namespace pde
} // namespace walberla
