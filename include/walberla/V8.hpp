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
//! \file V8.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

/**
 * @defgroup v8core V8 Core Library [Experimental]
 * @brief The waLBerla v8 core library.
 * 
 * The v8 core library is the new-generation portable framework core for waLBerla.
 * 
 * @warning This module is *under active development* and *highly experimental*.
 *          APIs may change rapidly and without warning.
 *          Use with care.
 * 
 * ## Include
 * 
 * ```
 * #include "walberla/V8.hpp"
 * ```
 */

#include "./v8/Memory.hpp"
#include "./v8/Sweep.hpp"
#include "./v8/HaloExchange.hpp"
#include "./v8/StencilRanges.hpp"
#include "./v8/Device.hpp"
