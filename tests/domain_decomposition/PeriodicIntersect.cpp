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
//! \file DataTypesTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/all.h"
#include "domain_decomposition/all.h"
#include "blockforest/all.h"

#include "core/math/AABB.h"
#include "core/math/Random.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"

#include "stencil/D3Q27.h"


namespace walberla {

bool periodicCheck( StructuredBlockForest& forest, const math::AABB& box1, const math::AABB& box2)
{
   bool flag = false;
   for (auto dir1 = stencil::D3Q27::begin(); dir1 != stencil::D3Q27::end(); ++dir1)
   {
      math::AABB tmp1 = box1;
      tmp1.translate( math::Vector3<real_t> (real_c(dir1.cx()) * forest.dx(), real_c(dir1.cy()) * forest.dy(), real_c(dir1.cz()) * forest.dz()));
      for (auto dir2 = stencil::D3Q27::begin(); dir2 != stencil::D3Q27::end(); ++dir2)
      {
         math::AABB tmp2 = box2;
         tmp2.translate( math::Vector3<real_t> (real_c(dir2.cx()) * forest.dx(), real_c(dir2.cy()) * forest.dy(), real_c(dir2.cz()) * forest.dz()));
         flag |= tmp1.intersects(tmp2);
         if (flag) break;
      }
      if (flag) break;
   }
   WALBERLA_CHECK_EQUAL( flag, forest.periodicIntersect(box1, box2), box1 << "\n" << box2 );
   return flag;
}

bool periodicCheck( StructuredBlockForest& forest, const math::AABB& box1, const math::AABB& box2, const real_t dx)
{
   bool flag = false;
   for (auto dir1 = stencil::D3Q27::begin(); dir1 != stencil::D3Q27::end(); ++dir1)
   {
      math::AABB tmp1 = box1;
      tmp1.translate( math::Vector3<real_t> (real_c(dir1.cx()) * forest.dx(), real_c(dir1.cy()) * forest.dy(), real_c(dir1.cz()) * forest.dz()));
      for (auto dir2 = stencil::D3Q27::begin(); dir2 != stencil::D3Q27::end(); ++dir2)
      {
         math::AABB tmp2 = box2;
         tmp2.translate( math::Vector3<real_t> (real_c(dir2.cx()) * forest.dx(), real_c(dir2.cy()) * forest.dy(), real_c(dir2.cz()) * forest.dz()));
         flag |= tmp1.intersects(tmp2, dx);
         if (flag) break;
      }
      if (flag) break;
   }
   WALBERLA_CHECK_EQUAL( flag, forest.periodicIntersect(box1, box2, dx), box1 << "\n" << box2 );
   return flag;
}

int main( int argc, char** argv )
{
   debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   WALBERLA_UNUSED(walberlaEnv);

   math::seedRandomGenerator(42);

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            false,                              // one block per process? - "false" means all blocks to one process
            true, true, true );                 // total periodicity

   math::AABB box1;
   math::AABB box2;

   // no periodicity check
   box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(0, 0, 0), math::Vector3<real_t>(1, 1, 1));
   real_t dx = real_c(0.1);
   for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
   {
      real_t delta = 0;

      delta = real_t(0.9);
      box2 = box1;
      box2.translate( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( box1.intersects( box2 ), box1 << "\n" << box2 );
      WALBERLA_CHECK( box1.intersects( box2, dx ), box1 << "\n" << box2 );

      delta = real_t(1.05);
      box2 = box1;
      box2.translate( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( !box1.intersects( box2 ), box1 << "\n" << box2 );
      WALBERLA_CHECK( box1.intersects( box2, dx ), box1 << "\n" << box2 );

      delta = real_t(1.15);
      box2 = box1;
      box2.translate( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( !box1.intersects( box2 ), box1 << "\n" << box2 );
      WALBERLA_CHECK( !box1.intersects( box2, dx ), box1 << "\n" << box2 );
   }

   // periodicity check
   auto temp = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(4, 4, 4), math::Vector3<real_t>(6, 6, 6));
   dx = real_c(0.3);
   for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
   {
      real_t delta = 0;

      delta = real_t(3);
      box1 = temp.getTranslated( -math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      box2 = temp.getTranslated( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( !forest->periodicIntersect(box1, box2), box1 << "\n" << box2 );
      WALBERLA_CHECK( !forest->periodicIntersect(box1, box2, dx), box1 << "\n" << box2 );

      delta = real_t(3.9);
      box1 = temp.getTranslated( -math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      box2 = temp.getTranslated( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( !forest->periodicIntersect(box1, box2), box1 << "\n" << box2 );
      WALBERLA_CHECK( forest->periodicIntersect(box1, box2, dx), box1 << "\n" << box2 );

      delta = real_t(4.1);
      box1 = temp.getTranslated( -math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      box2 = temp.getTranslated( math::Vector3<real_t> (real_c(dir.cx()) * delta, real_c(dir.cy()) * delta, real_c(dir.cz()) * delta));
      WALBERLA_CHECK( forest->periodicIntersect(box1, box2), box1 << "\n" << box2 );
      WALBERLA_CHECK( forest->periodicIntersect(box1, box2, dx), box1 << "\n" << box2 );
   }

   box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(1, 1, 1), math::Vector3<real_t>(2, 2, 2));
   box2 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(0, 0, 0), math::Vector3<real_t>(3, 3, 3));
   WALBERLA_CHECK( box1.intersects(box2) );
   WALBERLA_CHECK( forest->periodicIntersect(box1, box2) );
   periodicCheck(*forest, box1, box2);

   box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(1, 1, 1), math::Vector3<real_t>(2, 2, 2));
   box2 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(3, 3, 3), math::Vector3<real_t>(4, 4, 4));
   WALBERLA_CHECK( !box1.intersects(box2) );
   WALBERLA_CHECK( !forest->periodicIntersect(box1, box2) );
   periodicCheck(*forest, box1, box2);

   box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(0, 0, 0), math::Vector3<real_t>(2, 2, 2));
   box2 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(1, 1, 1), math::Vector3<real_t>(3, 3, 3));
   WALBERLA_CHECK( box1.intersects(box2) );
   WALBERLA_CHECK( forest->periodicIntersect(box1, box2) );
   periodicCheck(*forest, box1, box2);

   box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(-2, -2, -2), math::Vector3<real_t>(2, 2, 2));
   box2 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(9, 9, 9), math::Vector3<real_t>(11, 11, 11));
   WALBERLA_CHECK( !box1.intersects(box2) );
   WALBERLA_CHECK( forest->periodicIntersect(box1, box2) );
   periodicCheck(*forest, box1, box2);

   const int numTests = 100000;
   int count = 0;
   int countdx = 0;
   for (int i = 0; i < numTests; ++i)
   {
      box1 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(-10, -10, -10), math::Vector3<real_t>(math::realRandom<real_t>(-9,-5), math::realRandom<real_t>(-9,-5), math::realRandom<real_t>(-9,-5)));
      box2 = math::AABB::createFromMinMaxCorner( math::Vector3<real_t>(-10, -10, -10), math::Vector3<real_t>(math::realRandom<real_t>(-9,-5), math::realRandom<real_t>(-9,-5), math::realRandom<real_t>(-9,-5)));
      box1.translate( math::Vector3<real_t>( 2 * math::realRandom<real_t>() * forest->dx(), 2 * math::realRandom<real_t>() * forest->dy(), 2 * math::realRandom<real_t>() * forest->dz() ) );
      box2.translate( math::Vector3<real_t>( 2 * math::realRandom<real_t>() * forest->dx(), 2 * math::realRandom<real_t>() * forest->dy(), 2 * math::realRandom<real_t>() * forest->dz() ) );
      bool checkWithDx = periodicCheck( *forest, box1, box2, real_c(0.1) );
      if (checkWithDx)
      {
         ++countdx;
         WALBERLA_LOG_DETAIL("Intersection dx: " << countdx << " Periodic: " << !box1.intersects(box2));
         WALBERLA_UNUSED( countdx );
      }
      if (periodicCheck(*forest, box1, box2))
      {
         WALBERLA_CHECK(checkWithDx);
         ++count;
         WALBERLA_LOG_DETAIL("Intersection: " << count << " Periodic: " << !box1.intersects(box2));
         WALBERLA_UNUSED( count );
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}