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
//! \file Multigrid.impl.h
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "Multigrid.h"



namespace walberla {
namespace pde {



template< typename Stencil_T >
void Restrict< Stencil_T >::operator()( IBlock * const block ) const
{
   auto fine   = block->getData< Field_T >( fineFieldId_ );
   auto coarse = block->getData< Field_T >( coarseFieldId_ );

   WALBERLA_FOR_ALL_CELLS_XYZ( coarse,

      const cell_idx_t fx = 2*x;
      const cell_idx_t fy = 2*y;
      cell_idx_t fz = z;

      if( Stencil_T::D == uint_t(3) )
      {
         fz = 2*z;

         coarse->get(x,y,z) =   fine->get(fx  , fy  , fz+1) + fine->get(fx+1, fy  , fz+1)
                              + fine->get(fx  , fy+1, fz+1) + fine->get(fx+1, fy+1, fz+1);
      }
      else
      {
         WALBERLA_ASSERT_EQUAL(z, 0);
      }

      coarse->get(x,y,z) +=   fine->get(fx  , fy  , fz  ) + fine->get(fx+1, fy  , fz  )
                            + fine->get(fx  , fy+1, fz  ) + fine->get(fx+1, fy+1, fz  );
   );
}



template< typename Stencil_T >
void ProlongateAndCorrect< Stencil_T >::operator()( IBlock * const block ) const
{
   auto fine   = block->getData< Field_T >( fineFieldId_ );
   auto coarse = block->getData< Field_T >( coarseFieldId_ );

   WALBERLA_FOR_ALL_CELLS_XYZ( fine,
      if( Stencil_T::D == uint_t(3) )
      {
         fine->get(x,y,z) +=  real_t(0.125) * coarse->get(x/2,y/2,z/2);
      }
      else
      {
         WALBERLA_ASSERT_EQUAL(z, 0);
         fine->get(x,y,z) +=  real_t(0.25) * coarse->get(x/2,y/2,z);
      }
   );
}



template< typename Stencil_T >
void ComputeResidual< Stencil_T >::operator()( IBlock * const block ) const
{
   Field_T * rf = block->template getData< Field_T >( rId_            );
   Field_T * ff = block->template getData< Field_T >( fId_            );
   Field_T * uf = block->template getData< Field_T >( uId_            );
   StencilField_T * stencil = block->template getData< StencilField_T >( stencilId_ );

   WALBERLA_ASSERT_NOT_NULLPTR( rf      );
   WALBERLA_ASSERT_NOT_NULLPTR( ff      );
   WALBERLA_ASSERT_NOT_NULLPTR( uf      );
   WALBERLA_ASSERT_NOT_NULLPTR( stencil );

   WALBERLA_ASSERT_EQUAL( rf->xyzSize(), ff->xyzSize()      );
   WALBERLA_ASSERT_EQUAL( rf->xyzSize(), uf->xyzSize()      );
   WALBERLA_ASSERT_EQUAL( rf->xyzSize(), stencil->xyzSize() );

   WALBERLA_ASSERT_GREATER_EQUAL( uf->nrOfGhostLayers(), 1 );

   WALBERLA_FOR_ALL_CELLS_XYZ( uf,

      rf->get(x,y,z) = ff->get(x,y,z);

      for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
         rf->get(x,y,z) -= stencil->get( x, y, z, dir.toIdx() ) * uf->getNeighbor( x, y, z, *dir );
   )
}



template< typename Stencil_T >
void ComputeResidualFixedStencil< Stencil_T >::operator()( IBlock * const block ) const
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
         rf->get(x,y,z) -= weights_[ dir.toIdx() ] * uf->getNeighbor( x, y, z, *dir );
   )
}



} // namespace pde
} // namespace walberla
