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
//! \file ExtrapolationDirectionFinder.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#include "ExtrapolationDirectionFinder.h"

namespace walberla {
namespace pe_coupling {


void SphereNormalExtrapolationDirectionFinder
::getDirection( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const block, Vector3<cell_idx_t> & extrapolationDirection ) const
{
   BodyField_T * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );
   WALBERLA_ASSERT_NOT_NULLPTR( (*bodyField)(x,y,z) );

   real_t cx;
   real_t cy;
   real_t cz;
   blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

   Vector3<real_t> bodyCenterPosition = (*bodyField)(x,y,z)->getPosition();
   WALBERLA_ASSERT( !math::isnan(bodyCenterPosition) );

   Vector3<real_t> bodyNormal( cx - bodyCenterPosition[0], cy - bodyCenterPosition[1], cz - bodyCenterPosition[2] );

   findCorrespondingLatticeDirection< stencil::D3Q27 >( bodyNormal, extrapolationDirection );
}

} // namespace pe_coupling
} // namespace walberla
