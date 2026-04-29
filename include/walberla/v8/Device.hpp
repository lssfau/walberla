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
//! \file Device.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

/**
 * @weakgroup v8core-device Device-Accessible APIs
 * @brief Classes and functions available to device kernels
 * @ingroup v8core
 * 
 * This category lists classes and functions from both the `walberla::v8` modules
 * and the legacy `core` module that are available to CUDA/HIP device code.
 * 
 * @warning These APIs are only available to device code if the V8 core library is enabled in your build.
 */

#include "core/Stdlib.hpp"
