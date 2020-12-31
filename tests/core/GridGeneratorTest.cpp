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
//! \file GridGeneratorTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include "core/grid_generator/HCPIterator.h"
#include "core/grid_generator/SCIterator.h"

using namespace walberla;
using namespace walberla::grid_generator;

void unitCellTest()
{
   using namespace walberla::debug::check_functions_detail;

   bool checkpoint, origin;
   const real_t spacing = real_c(1.0);
   const math::AABB aabb(-5 * spacing, -5 * spacing, -5 * spacing, 5 * spacing, 5 * spacing, 5 * spacing);

   checkpoint = false;
   origin = false;
   for (auto it = HCPIterator(aabb, Vector3<real_t>(-HCPIterator::getUnitCellX(spacing),0,0), spacing); it != HCPIterator(); ++it)
   {
      if (check_float_equal(*it, Vector3<real_t>(0,0,0))) checkpoint = true;
      if (check_float_equal(*it, Vector3<real_t>(-HCPIterator::getUnitCellX(spacing),0,0))) origin = true;
   }
   WALBERLA_CHECK(checkpoint);
   WALBERLA_CHECK(origin);

   checkpoint = false;
   origin = false;
   for (auto it = HCPIterator(aabb, Vector3<real_t>(0,-HCPIterator::getUnitCellY(spacing),0), spacing); it != HCPIterator(); ++it)
   {
      if (check_float_equal(*it, Vector3<real_t>(0,0,0))) checkpoint = true;
      if (check_float_equal(*it, Vector3<real_t>(0,-HCPIterator::getUnitCellY(spacing),0))) origin = true;
   }
   WALBERLA_CHECK(checkpoint);
   WALBERLA_CHECK(origin);

   checkpoint = false;
   origin = false;
   for (auto it = HCPIterator(aabb, Vector3<real_t>(0,0,-HCPIterator::getUnitCellZ(spacing)), spacing); it != HCPIterator(); ++it)
   {
      if (check_float_equal(*it, Vector3<real_t>(0,0,0))) checkpoint = true;
      if (check_float_equal(*it, Vector3<real_t>(0,0,-HCPIterator::getUnitCellZ(spacing)))) origin = true;
   }
   WALBERLA_CHECK(checkpoint);
   WALBERLA_CHECK(origin);
}

template <class GridGenerator>
void referencePointTest()
{
   math::AABB domain(0,0,0,10,10,10);
   real_t spacing = real_c(1);
   auto dx = Vector3<real_t>(GridGenerator::getUnitCellX(spacing), GridGenerator::getUnitCellY(spacing), GridGenerator::getUnitCellZ(spacing) );
   auto lowerIt  = GridGenerator(domain, Vector3<real_t>(5,5,5) - dx * 30, spacing);
   auto inIt     = GridGenerator(domain, Vector3<real_t>(5,5,5), spacing);
   auto higherIt = GridGenerator(domain, Vector3<real_t>(5,5,5) + dx * 30, spacing);
   auto endIt = GridGenerator();
   for ( ; lowerIt != endIt; ++lowerIt, ++inIt, ++higherIt)
   {
      WALBERLA_CHECK_FLOAT_EQUAL( *lowerIt, *inIt);
      WALBERLA_CHECK_FLOAT_EQUAL( *higherIt, *inIt);
   }
   WALBERLA_CHECK_EQUAL( lowerIt, endIt);
   WALBERLA_CHECK_EQUAL( inIt, endIt);
   WALBERLA_CHECK_EQUAL( higherIt, endIt);
}

template <class Grid>
void rangeBasedTest()
{
   math::AABB domain(0,0,0,10,10,10);
   real_t spacing = real_c(1);
   auto dx = Vector3<real_t>( Grid::iterator::getUnitCellX(spacing), Grid::iterator::getUnitCellY(spacing), Grid::iterator::getUnitCellZ(spacing) );
   auto lowerIt  = typename Grid::iterator(domain, Vector3<real_t>(5,5,5) - dx * 30, spacing);
   auto endIt = typename Grid::iterator();
   for ( const auto pt : Grid(domain, Vector3<real_t>(5,5,5) - dx * 30, spacing) )
   {
      WALBERLA_CHECK( lowerIt != endIt );
      WALBERLA_CHECK_FLOAT_EQUAL( *lowerIt, pt);
      ++lowerIt;
   }
   WALBERLA_CHECK( lowerIt == endIt );
}

int main( int argc, char** argv )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   std::vector< Vector3<real_t> > points;
   points.emplace_back( real_t(0), real_t(0), real_t(0) );
   points.emplace_back( real_t(1), real_t(0), real_t(0) );
   points.emplace_back( real_t(0), real_t(1), real_t(0) );
   points.emplace_back( real_t(1), real_t(1), real_t(0) );
   points.emplace_back( real_t(0), real_t(0), real_t(1) );
   points.emplace_back( real_t(1), real_t(0), real_t(1) );
   points.emplace_back( real_t(0), real_t(1), real_t(1) );
   points.emplace_back( real_t(1), real_t(1), real_t(1) );
   auto correctPointIt = points.begin();
   for (auto it = SCIterator(AABB(real_c(-0.01), real_c(-0.01), real_c(-0.01), real_c(1.9),real_c(1.9),real_c(1.9)), Vector3<real_t>(0,0,0), 1); it != SCIterator(); ++it, ++correctPointIt)
      WALBERLA_CHECK_FLOAT_EQUAL( *it, *correctPointIt, (*it) << "!=" << (*correctPointIt) );

   points.clear();
   points.emplace_back(real_t(0),real_t(0),real_t(0) );
   points.emplace_back(real_t(1),real_t(0),real_t(0) );
   points.emplace_back(real_c(0.5), real_c(0.866025),real_t(0) );
   points.emplace_back(real_c(1.5), real_c(0.866025),real_t(0) );
   points.emplace_back(real_t(0), real_c(1.73205),real_t(0) );
   points.emplace_back(real_t(1), real_c(1.73205),real_t(0) );

   points.emplace_back(real_c(0.5), real_c(0.288675), real_c(0.816497) );
   points.emplace_back(real_c(1.5), real_c(0.288675), real_c(0.816497) );
   points.emplace_back(real_t(0), real_c(1.1547), real_c(0.816497) );
   points.emplace_back(real_t(1), real_c(1.1547), real_c(0.816497) );

   points.emplace_back(real_t(0),real_t(0), real_c(1.63299) );
   points.emplace_back(real_t(1),real_t(0), real_c(1.63299) );
   points.emplace_back(real_c(0.5), real_c(0.866025), real_c(1.63299) );
   points.emplace_back(real_c(1.5), real_c(0.866025), real_c(1.63299) );
   points.emplace_back(real_t(0), real_c(1.73205), real_c(1.63299) );
   points.emplace_back(real_t(1), real_c(1.73205), real_c(1.63299) );
   correctPointIt = points.begin();
   for (auto it = HCPIterator(AABB(real_c(-0.01), real_c(-0.01), real_c(-0.01), real_c(1.9),real_c(1.9),real_c(1.9)), Vector3<real_t>(real_t(0),real_t(0),real_t(0)), 1); it != HCPIterator(); ++it, ++correctPointIt)
   {
      WALBERLA_CHECK( floatIsEqual((*it)[0], (*correctPointIt)[0], real_c(0.00001)), (*it) << "!=" << (*correctPointIt));
      WALBERLA_CHECK( floatIsEqual((*it)[1], (*correctPointIt)[1], real_c(0.00001)), (*it) << "!=" << (*correctPointIt));
      WALBERLA_CHECK( floatIsEqual((*it)[2], (*correctPointIt)[2], real_c(0.00001)), (*it) << "!=" << (*correctPointIt));
   }

   unitCellTest();
   referencePointTest<SCIterator>();
   referencePointTest<HCPIterator>();

   rangeBasedTest<SCGrid>();
   rangeBasedTest<HCPGrid>();

   return EXIT_SUCCESS;
}
