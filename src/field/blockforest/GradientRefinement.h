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
//! \file GradientRefinement.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockForest.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "field/iterators/IteratorMacros.h"



namespace walberla {
namespace field {



template< typename VectorField_T, typename Filter_T, bool Pseudo2D = false >
class GradientRefinement
{
public:

   GradientRefinement( const ConstBlockDataID & fieldId, const Filter_T & filter,
                       const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel ) :
      fieldId_( fieldId ), filter_( filter ),
      upperLimit_( upperLimit ), lowerLimit_( lowerLimit ), maxLevel_( maxLevel )
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const BlockForest & forest );

private:

   void calculate( const Vector3<real_t> & v, const Vector3<real_t> & vn, real_t * gradient )
   {
      gradient[0] = std::abs( v[0] - vn[0] );
      gradient[1] = std::abs( v[1] - vn[1] );
      gradient[2] = std::abs( v[2] - vn[2] );
   }

   ConstBlockDataID fieldId_;

   Filter_T filter_;

   real_t upperLimit_;
   real_t lowerLimit_;

   uint_t maxLevel_;

}; // class GradientRefinement



template< typename VectorField_T, typename Filter_T, bool Pseudo2D >
void GradientRefinement< VectorField_T, Filter_T, Pseudo2D >::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                                                          std::vector< const Block * > &, const BlockForest & )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const Block * const block = it->first;
      const VectorField_T * u = block->template getData< VectorField_T >( fieldId_ );

      if( u == nullptr )
      {
         it->second = uint_t(0);
         continue;
      }

      CellInterval interval = u->xyzSize();
      CellInterval innerInterval = interval;
      Cell expand( cell_idx_t(-1), cell_idx_t(-1), Pseudo2D ? cell_idx_t(0) : cell_idx_t(-1) );
      innerInterval.expand( expand );

      const int gradients = Pseudo2D ? 12 : 18;

      const cell_idx_t one( cell_idx_t(1) );

      bool refine( false );
      bool coarsen( true );

      filter_( *block );

      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,

         if( filter_(x,y,z) && !refine )
         {
            const Vector3< real_t > & v = u->get(x,y,z);

            real_t gradient[gradients];
            for( int i = 0; i < gradients; ++i )
               gradient[i] = real_t(0);

            if( innerInterval.contains(x,y,z) )
            {
               if( filter_(x+one,y,z) )
               {
                  const Vector3< real_t > & vn = u->get(x+one,y,z);
                  calculate( v, vn, gradient + 0 );
               }
               if( filter_(x-one,y,z) )
               {
                  const Vector3< real_t > & vn = u->get(x-one,y,z);
                  calculate( v, vn, gradient + 3 );
               }
               if( filter_(x,y+one,z) )
               {
                  const Vector3< real_t > & vn = u->get(x,y+one,z);
                  calculate( v, vn, gradient + 6 );
               }
               if( filter_(x,y-one,z) )
               {
                  const Vector3< real_t > & vn = u->get(x,y-one,z);
                  calculate( v, vn, gradient + 9 );
               }
               if( !Pseudo2D )
               {
                  if( filter_(x,y,z+one) )
                  {
                     const Vector3< real_t > & vn = u->get(x,y,z+one);
                     calculate( v, vn, gradient + 12 );
                  }
                  if( filter_(x,y,z-one) )
                  {
                     const Vector3< real_t > & vn = u->get(x,y,z-one);
                     calculate( v, vn, gradient + 15 );
                  }
               }
            }
            else
            {
               if( interval.contains(x+one,y,z) && filter_(x+one,y,z) )
               {
                  const Vector3< real_t > & vn = u->get(x+one,y,z);
                  calculate( v, vn, gradient + 0 );
               }
               if( interval.contains(x-one,y,z) && filter_(x-one,y,z) )
               {
                  const Vector3< real_t > & vn = u->get(x-one,y,z);
                  calculate( v, vn, gradient + 3 );
               }
               if( interval.contains(x,y+one,z) && filter_(x,y+one,z) )
               {
                  const Vector3< real_t > & vn = u->get(x,y+one,z);
                  calculate( v, vn, gradient + 6 );
               }
               if( interval.contains(x,y-one,z) && filter_(x,y-one,z) )
               {
                  const Vector3< real_t > & vn = u->get(x,y-one,z);
                  calculate( v, vn, gradient + 9 );
               }
               if( !Pseudo2D )
               {
                  if( interval.contains(x,y,z+one) && filter_(x,y,z+one) )
                  {
                     const Vector3< real_t > & vn = u->get(x,y,z+one);
                     calculate( v, vn, gradient + 12 );
                  }
                  if( interval.contains(x,y,z-one) && filter_(x,y,z-one) )
                  {
                     const Vector3< real_t > & vn = u->get(x,y,z-one);
                     calculate( v, vn, gradient + 15 );
                  }
               }
            }

            real_t magnitute( real_t(0) );
            for( int i = 0; i < gradients; ++i )
               magnitute += gradient[i];

            if( magnitute > lowerLimit_ )
            {
               coarsen = false;
               if( magnitute > upperLimit_ )
                  refine = true;
            }
         }
      )

      if( refine && block->getLevel() < maxLevel_ )
      {
         WALBERLA_ASSERT( !coarsen );
         it->second = block->getLevel() + uint_t(1);
      }
      if( coarsen && block->getLevel() > uint_t(0) )
      {
         WALBERLA_ASSERT( !refine );
         it->second = block->getLevel() - uint_t(1);
      }
   }
}



} // namespace field
} // namespace walberla
