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
//! \file HaloExchange.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

/**
 * @defgroup v8core-haloexchange Halo Exchange
 * @ingroup v8core
 * @brief Classes facilitating exchange of halo-layer data between blocks.
 * 
 * ## Include
 * 
 * ```
 * #include "walberla/v8/HaloExchange.hpp"
 * ```
 */

#include "./halo_exchange/HaloExchange.hpp"
#include "./halo_exchange/LegacyCpuScheme.hpp"
#include "./halo_exchange/LegacyGpuScheme.hpp"

#include "./halo_exchange/fields/GenericFieldPackInfos.hpp"
#include "./halo_exchange/fields/PackInfoSelection.hpp"
