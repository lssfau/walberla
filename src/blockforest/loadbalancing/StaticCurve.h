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
//! \file StaticCurve.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/SetupBlockForest.h"



namespace walberla {
namespace blockforest {



class StaticLevelwiseCurveBalance // all blocks are assumed to have the same weight/workload
{
public:

   StaticLevelwiseCurveBalance( const bool hilbert = true ) : hilbert_( hilbert ) {}

   uint_t operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t perProcessMemoryLimit );

private:

   bool hilbert_;
};



class StaticLevelwiseCurveBalanceWeighted // takes the weight/workload of all blocks into account
{
public:

   StaticLevelwiseCurveBalanceWeighted( const bool hilbert = true ) : hilbert_( hilbert ) {}

   uint_t operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t perProcessMemoryLimit );

private:

   bool hilbert_;
};



} // namespace blockforest
} // namespace walberla
