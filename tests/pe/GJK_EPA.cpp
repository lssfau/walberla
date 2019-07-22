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
//! \file GJK_EPA.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/all.h"
#include "core/Abort.h"
#include "core/timing/Timer.h"
#include "core/math/Sample.h"
#include "core/mpi/Environment.h"
#include "stencil/D3Q27.h"
#include "pe/Types.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/collision/Collide.h"

#include "pe/collision/GJKEPAHelper.h"

#include <cstdlib>

using namespace walberla;
using namespace walberla::pe;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   real_t radius            = real_c(1);
   real_t overlap           = real_c(0.2);
   unsigned long iterations = 100;
   real_t epaTolerance      = real_t(0.000001);

   if (argc != 5)
   {
      WALBERLA_ABORT("usage: ./PE_GJKEPA radius overlap maxIterations epaTolerance");
   }
   radius       = atof(argv[1]);
   overlap      = atof(argv[2]);
   iterations   = atol(argv[3]);
   epaTolerance = atof(argv[4]);

   WALBERLA_LOG_INFO("radius       : " << radius);
   WALBERLA_LOG_INFO("overlap      : " << overlap << " of radius");
   WALBERLA_LOG_INFO("gjkIterations: " << iterations);
   WALBERLA_LOG_INFO("epaTolerance : " << epaTolerance);


   Sphere s0( 0, 0, Vec3( 0, 0, 0 ), Vec3(), Quat(), real_t( radius ), 0, false, true, false );
   Sphere s1( 1, 0, Vec3( 1, 1, 1 ).getNormalized() * ( real_t(2) - overlap ) * radius, Vec3(), Quat(), real_t( radius ), 0, false, true, false );

   Vec3   contactPoint;
   Vec3   normal;
   real_t penetrationDepth;

   Vec3   contactPointGJK;
   Vec3   normalGJK;
   real_t penetrationDepthGJK;

   math::Sample point;
   math::Sample norm;
   math::Sample pen;

   int iter = 1000;
   WcTimer timer;
   timer.start();
   timer.end();
   while (timer.last() < real_t(1))
   {
      iter *= 2;
      timer.start();
      for (int i = 0; i < iter; ++i)
         collide(&s0, &s1, contactPoint, normal, penetrationDepth);
      timer.end();
   }
   WALBERLA_LOG_RESULT("analytical: " << timer.last() << "s for " << iter << " iterations");
   double timeAnalytical = timer.last() / iter;

   iter = 1000;
   timer.start();
   timer.end();
   while (timer.last() < real_t(1))
   {
      iter *= 2;
      timer.start();
      for (int i = 0; i < iter; ++i)
         collideGJK(&s0, &s1, contactPointGJK, normalGJK, penetrationDepthGJK, iterations, epaTolerance);
      timer.end();
   }
   WALBERLA_LOG_RESULT("GJK: " << timer.last() << "s for " << iter << " iterations");
   WALBERLA_LOG_RESULT("SLOWDOWN: " << timer.last() / iter / timeAnalytical << "x");


   for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
   {
      Vec3 pos = Vec3( real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()) ).getNormalized() * ( real_t(2) - overlap ) * radius;
      s1.setPosition(pos);

      if (!collide(&s0, &s1, contactPoint, normal, penetrationDepth))
      {
         WALBERLA_ABORT("no collision detected! (analytical)");
      }
      if (!collideGJK(&s0, &s1, contactPointGJK, normalGJK, penetrationDepthGJK, iterations, epaTolerance))
      {
         WALBERLA_LOG_WARNING("position        : " << pos);
         WALBERLA_LOG_WARNING("contactPoint    : " << contactPoint);
         WALBERLA_LOG_WARNING("contactNormal   : " << normal);
         WALBERLA_LOG_WARNING("penetrationDepth: " << penetrationDepth);
         WALBERLA_ABORT("no collision detected! (GJK & EPA)");
      }

      WALBERLA_ASSERT_FLOAT_EQUAL( normalGJK.length(), real_t(1));

//      WALBERLA_LOG_INFO("**** (" << dir.cx() << ", " << dir.cy() << ", " << dir.cz() << ") ****");
//      WALBERLA_LOG_INFO("deltaPoint : |" << contactPoint - contactPointGJK << "| = " << (contactPoint - contactPointGJK).length() );
//      WALBERLA_LOG_INFO("deltaNormal: |" << normal - normalGJK << "| = " << (normal - normalGJK).length() << " (" << acos(normal*normalGJK) / math::pi * 180 << "°)" );
//      WALBERLA_LOG_INFO("deltaPen   : " << penetrationDepth << " - " << penetrationDepthGJK << " = " << penetrationDepth - penetrationDepthGJK);

      point.insert( (contactPoint - contactPointGJK).length() );
      norm.insert( acos(normal*normalGJK) / math::pi * 180 );
      pen.insert( fabs(penetrationDepth - penetrationDepthGJK) );
   }

   WALBERLA_LOG_RESULT( "max contact point deviation    : " << point.max() );
   WALBERLA_LOG_RESULT( "max contact normal deviation   : " << norm.max() << "°" );
   WALBERLA_LOG_RESULT( "max penetration depth deviation: " << pen.max() );

   return EXIT_SUCCESS;
}
