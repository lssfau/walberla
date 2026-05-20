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
//! \file CommSchemeOptions.hpp
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Set.h"
#include "core/uid/SUID.h"

namespace walberla::v8::halo_exchange
{

/**
 * @brief Options controlling construction of the underlying comm scheme.
 * @ingroup v8core-haloexchange
 * 
 * Options that are not applicable to the selected comm scheme are ignored.
 * 
 * Supported comm schemes:
 *  - `UniformBufferedScheme`
 *  - `UniformGPUScheme`
 * 
 */
struct CommSchemeOptions
{
   Set< SUID > requiredBlockSelectors     = Set< SUID >::emptySet();
   Set< SUID > incompatibleBlockSelectors = Set< SUID >::emptySet();
   bool sendDirectlyFromGPU               = false;
   bool useLocalCommunication             = true;
   std::optional< int > mpiTag;
};

} // namespace walberla::v8::halo_exchange
