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
//! \file Sanitizer.h
//! \ingroup core
//! \author Dominik Thoennes <dominik.thoennes@fau.de>
//
//======================================================================================================================

#pragma once

#if defined(WALBERLA_CXX_COMPILER_IS_CLANG) || defined(WALBERLA_CXX_COMPILER_IS_GNU)
# define ATTRIBUTE_NO_SANITIZE_ADDRESS __attribute__((no_sanitize_address))
#else
# define ATTRIBUTE_NO_SANITIZE_ADDRESS
#endif

#if defined(WALBERLA_CXX_COMPILER_IS_GNU)
# define ATTRIBUTE_NO_SANITIZE_UNDEFINED __attribute__((no_sanitize_undefined))
#elif defined(WALBERLA_CXX_COMPILER_IS_CLANG)
# define ATTRIBUTE_NO_SANITIZE_UNDEFINED __attribute__((no_sanitize("undefined")))
#else
# define ATTRIBUTE_NO_SANITIZE_UNDEFINED
#endif