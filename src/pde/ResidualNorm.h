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
//! \file ResidualNorm.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/Reduce.h"
#include "domain_decomposition/BlockStorage.h"
#include "field/GhostLayerField.h"
#include "field/iterators/IteratorMacros.h"



namespace walberla {
namespace pde {



template< typename Stencil_T >
class ResidualNorm
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;
   
   ResidualNorm( const BlockStorage & blocks, const ConstBlockDataID & uId, const ConstBlockDataID & fId, const std::vector< real_t > & weights,
                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), uId_( uId ), fId_( fId ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {
      WALBERLA_ASSERT_EQUAL( weights.size(), Stencil_T::Size );
      for( uint_t i = uint_t(0); i < Stencil_T::Size; ++i )
         weights_[i] = weights[i];
      init();
   }

   real_t operator()() const { return weightedL2(); }

   real_t weightedL2() const;
   
protected:

   void init();

   

   const BlockStorage & blocks_;
   
   ConstBlockDataID uId_;
   ConstBlockDataID fId_;
   
   real_t weights_[ Stencil_T::Size ];
   
   real_t cells_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



template< typename Stencil_T >
real_t ResidualNorm< Stencil_T >::weightedL2() const
{
   real_t result( real_t(0) );
   
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      const Field_T * const uf = block->template getData< const Field_T >( uId_ );
      const Field_T * const ff = block->template getData< const Field_T >( fId_ );
      
      real_t blockResult( real_t(0) );
      
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP( uf, omp parallel for schedule(static) reduction(+:blockResult),

         real_t d = ff->get(x,y,z);

         for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
            d -= weights_[ dir.toIdx() ] * uf->getNeighbor(x,y,z,*dir);

         d -= weights_[ Stencil_T::idx[ stencil::C ] ] * uf->get(x,y,z);
         
         blockResult += d * d;

      ) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

      result += blockResult;
   }
   
   mpi::allReduceInplace( result, mpi::SUM );
   return std::sqrt( result / cells_ );
}



template< typename Stencil_T >
void ResidualNorm< Stencil_T >::init()
{
   uint_t cells( uint_t(0) );
   
   for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
   {
      const Field_T * const u = block->template getData< const Field_T >( uId_ );
      cells += u->xyzSize().numCells();
   }
   
   cells_ = real_c( cells );
   mpi::allReduceInplace( cells_, mpi::SUM );
}



} // namespace pde
} // namespace walberla
