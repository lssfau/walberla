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
//! \file HCSITS.cpp
//! \brief checks equality of hash grids and simple ccd
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"

#include "blockforest/Initialization.h"
#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   // create blocks
   auto forest = blockforest::createBlockForest( math::AABB(0,0,0,10,10,10),
                                                 Vector3<uint_t>(1,1,1),
                                                 Vector3<bool>(false, false, false),
                                                 1,
                                                 1);

   WALBERLA_CHECK_EQUAL( forest->size(), 8);
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      blockforest::Block& currentBlock = *(dynamic_cast<blockforest::Block*>(&(*blockIt)));
      WALBERLA_CHECK_EQUAL( currentBlock.getLevel(), 1);
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
