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
//! \file all.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Collective header file for module core
//
//======================================================================================================================

#pragma once

#include "Abort.h"
#include "AllSet.h"
#include "Array.h"
#include "DataTypes.h"
#include "Deprecated.h"
#include "EndianIndependentSerialization.h"
#include "Environment.h"
#include "FunctionTraits.h"
#include "GetPID.h"
#include "Hostname.h"
#include "Macros.h"
#include "NonCopyable.h"
#include "NonCreateable.h"
#include "OpenMP.h"
#include "Regex.h"
#include "Set.h"
#include "SharedFunctor.h"
#include "Sleep.h"
#include "StringUtility.h"
#include "VectorTrait.h"

#include "cell/all.h"
#include "config/all.h"
#include "debug/all.h"
#include "grid_generator/all.h"
#include "logging/all.h"
#include "math/all.h"
#include "mpi/all.h"
#include "perf_analysis/all.h"
#include "selectable/all.h"
#include "singleton/all.h"
#include "timing/all.h"
#include "uid/all.h"
