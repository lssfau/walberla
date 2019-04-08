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

#include "stencil/D3Q7.h"

#include <array>

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
         coarse->get(x,y,z) = real_t(0);
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



template< typename Stencil_T >
void CoarsenStencilFieldsDCA<Stencil_T >::operator()( const std::vector<BlockDataID> & stencilFieldId ) const
{
   const real_t scalingFactor = real_t(1)/real_c(2<< (operatorOrder_-1) ); // scaling by ( 1/h^operatorOrder )^lvl

   WALBERLA_ASSERT_EQUAL(numLvl_, stencilFieldId.size(), "This function can only be called when operating with stencil fields!");

   for ( uint_t lvl = 1; lvl < numLvl_; ++lvl )
   {
      for( auto block = blocks_->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_->end(); ++block )
      {
         StencilField_T * fine   = block->template getData< StencilField_T >( stencilFieldId[lvl-1] );
         StencilField_T * coarse = block->template getData< StencilField_T >( stencilFieldId[lvl] );

         WALBERLA_FOR_ALL_CELLS_XYZ(coarse,
            for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
               coarse->get(x,y,z, dir.toIdx()) = scalingFactor * fine->get(x,y,z, dir.toIdx());
         )
      }
   }
}



template< >
void CoarsenStencilFieldsGCA< stencil::D3Q7 >::operator()( const std::vector<BlockDataID> & stencilFieldId ) const
{

   WALBERLA_ASSERT_EQUAL(numLvl_, stencilFieldId.size(), "This function can only be called when operating with stencil fields!");

   // Apply Galerkin coarsening to each level //
   // currently only implemented for CCMG with constant restriction and prolongation

   for ( uint_t lvl = 1; lvl < numLvl_; ++lvl )
   {
      for( auto block = blocks_->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_->end(); ++block )
      {
         StencilField_T * fine   = block->getData< StencilField_T >( stencilFieldId[lvl-1] );
         StencilField_T * coarse = block->getData< StencilField_T >( stencilFieldId[lvl] );


         typedef std::array<std::array<std::array<real_t,7>,7>,7> Array3D_7;
         typedef std::array<std::array<std::array<real_t,2>,2>,2> Array3D_2;

         Array3D_2 r;
         Array3D_7 p;

         // Set to restriction weights from constant restrict operator
          for(auto z = uint_t(0); z < uint_t(2); ++z) {
            for(auto y = uint_t(0); y < uint_t(2); ++y) {
               for(auto x = uint_t(0); x < uint_t(2); ++x) {
                  r[x][y][z] = real_t(1);
               }
            }
         }

         // Init to 0 //
         for(auto k = uint_t(0); k < uint_t(7); ++k) {
            for(auto j = uint_t(0); j < uint_t(7); ++j) {
               for(auto i = uint_t(0); i < uint_t(7); ++i) {
                  p[i][j][k] = real_t(0);
               }
            }
         }

         // Set center to prolongation weights, including overrelaxation factor (latter therefore no longer needed in ProlongateAndCorrect)
         for(auto k = uint_t(0); k < uint_t(2); ++k) {
            for(auto j = uint_t(0); j < uint_t(2); ++j) {
               for(auto i = uint_t(0); i < uint_t(2); ++i) {
                  p[i+2][j+2][k+2] = real_t(0.125)/overrelaxFact_;   // Factor 0.125 such that same prolongation operator for DCA and GCA
               }
            }
         }


         WALBERLA_FOR_ALL_CELLS_XYZ(coarse,

            Array3D_7 ap;

            const cell_idx_t fx = 2*x;
            const cell_idx_t fy = 2*y;
            const cell_idx_t fz = 2*z;

            for(auto k = uint_t(0); k < uint_t(7); ++k)
               for(auto j = uint_t(0); j < uint_t(7); ++j)
                  for(auto i = uint_t(0); i < uint_t(7); ++i)
                     ap[i][j][k] = real_t(0);


            // Tested for spatially varying stencils! //
            for(auto k = uint_t(1); k < uint_t(5); ++k)
               for(auto j = uint_t(1); j < uint_t(5); ++j)
                  for(auto i = uint_t(1); i < uint_t(5); ++i){
                     ap[i][j][k] = real_t(0);
                     for(auto d = stencil::D3Q7::begin(); d != stencil::D3Q7::end(); ++d ){
                        ap[i][j][k] += p[ uint_c(cell_idx_c(i)+d.cx()) ] [ uint_c(cell_idx_c(j)+d.cy()) ] [uint_c(cell_idx_c(k)+d.cz()) ] * fine->get( fx+cell_idx_c(i%2), fy+cell_idx_c(j%2), fz+cell_idx_c(k%2), d.toIdx() ); // contains elements of one row of coarse grid matrix
                     }
                  }

            // Checked, correct! //
            for(auto d = stencil::D3Q7::begin(); d != stencil::D3Q7::end(); ++d ){
               real_t sum = 0;
               for(auto k = uint_t(0); k < uint_t(2); ++k)
                  for(auto j = uint_t(0); j < uint_t(2); ++j)
                     for(auto i = uint_t(0); i < uint_t(2); ++i) {
                        sum += ap[uint_c(cell_idx_c(i)+2-2*d.cx())] [uint_c(cell_idx_c(j)+2-2*d.cy())] [uint_c(cell_idx_c(k)+2-2*d.cz())] * r[i][j][k];
                        // either i+0 or i+4    // either j+0 or j+4    // either k+0 or k+4       // always 1 here
                     }

               coarse->get(x,y,z,*d) = sum;
            }
         )
      }
   }
}


} // namespace pde
} // namespace walberla
