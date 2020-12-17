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
//! \file python/Python.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

// Do not reorder includes - the include order is important
#include "waLBerlaDefinitions.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON // macro defined in waLBerlaDefinitions.h

#ifdef _MSC_VER
#pragma warning ( push, 3 )
#pragma warning ( disable: 4244 4275 4800 4251 4267 )
#ifndef HAVE_ROUND
#define HAVE_ROUND 1
#define __CREATED_HAVE_ROUND
#endif
#endif

#include "pybind11/pybind11.h"

#ifdef _MSC_VER
#ifdef __CREATED_HAVE_ROUND
#undef HAVE_ROUND
#undef __CREATED_HAVE_ROUND
#endif
#pragma warning ( pop )
#endif

#endif
