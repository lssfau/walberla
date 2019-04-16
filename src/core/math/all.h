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

#include "AABB.h"
#include "Constants.h"
#include "DistributedSample.h"
#include "FPClassify.h"
#include "FastInvSqrt.h"
#include "GenericAABB.h"
#include "IntegerFactorization.h"
#include "KahanSummation.h"
#include "Limits.h"
#include "MathTrait.h"
#include "Matrix2.h"
#include "Matrix3.h"
#include "Parser.h"
#include "ParserOMP.h"
#ifdef WALBERLA_BUILD_WITH_BOOST
#include "PhysicalCheck.h"
#endif
#include "Plane.h"
#include "Primes.h"
#include "Random.h"
#include "Sample.h"
#include "SqrtTrait.h"
#include "Uint.h"
#include "Utility.h"
#include "Vector2.h"
#include "Vector3.h"

#include "equation_system/all.h"
