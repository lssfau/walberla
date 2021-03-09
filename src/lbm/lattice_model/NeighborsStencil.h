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
//! \file NeighborsStencil.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"

#include <type_traits>



namespace walberla {
namespace lbm {



template< typename LatticeModel_T, class Enable = void  >
struct NeighborsStencil;

template< typename LatticeModel_T >
struct NeighborsStencil< LatticeModel_T, typename std::enable_if< LatticeModel_T::Stencil::D == 2 >::type >
{
   using type = stencil::D2Q9;
};

template< typename LatticeModel_T >
struct NeighborsStencil< LatticeModel_T, typename std::enable_if< LatticeModel_T::Stencil::D == 3 >::type >
{
   using type = stencil::D3Q27;
};



} // namespace lbm
} // namespace walberla
