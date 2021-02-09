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
//! \file Cartesian.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/SetupBlockForest.h"



namespace walberla {
namespace blockforest {



class CartesianDistribution
{
public:

   CartesianDistribution( const uint_t numberOfXProcesses, const uint_t numberOfYProcesses, const uint_t numberOfZProcesses,
                          std::vector< uint_t > * processIdMap = nullptr ) :
      numberOfXProcesses_( numberOfXProcesses ), numberOfYProcesses_( numberOfYProcesses ), numberOfZProcesses_( numberOfZProcesses ),
      processIdMap_( processIdMap )
   {}

   uint_t operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t perProcessMemoryLimit );

private:

   uint_t numberOfXProcesses_;
   uint_t numberOfYProcesses_;
   uint_t numberOfZProcesses_;

   std::vector< uint_t > * processIdMap_;
};



} // namespace blockforest
} // namespace walberla
