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
//! \file MapToPeriodicDomain.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/MapPointToPeriodicDomain.h"

namespace walberla {
namespace testing {

int main (int argc, char** argv)
{
   using namespace walberla::domain_decomposition;

   debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   WALBERLA_UNUSED(walberlaEnv);

   std::array<bool, 3> isPeriodic{{true, true, true}};
   math::Vector3<real_t> minCorner(real_t(2), real_t(3), real_t(4));
   math::Vector3<real_t> maxCorner(real_t(3), real_t(4), real_t(5));
   math::AABB domain(minCorner, maxCorner);

   //minCorner -> minCorner
   math::Vector3<real_t> p = minCorner;
   mapPointToPeriodicDomain( isPeriodic, domain, p);
   WALBERLA_CHECK_IDENTICAL( p[0], minCorner[0] );
   WALBERLA_CHECK_IDENTICAL( p[1], minCorner[1] );
   WALBERLA_CHECK_IDENTICAL( p[2], minCorner[2] );

   //maxCorner -> minCorner
   p = maxCorner;
   mapPointToPeriodicDomain( isPeriodic, domain, p);
   WALBERLA_CHECK_IDENTICAL( p[0], minCorner[0] );
   WALBERLA_CHECK_IDENTICAL( p[1], minCorner[1] );
   WALBERLA_CHECK_IDENTICAL( p[2], minCorner[2] );

   //minCorner - epsilon -> maxCorner - epsilon
   p[0] = std::nextafter(domain.xMin(), real_t(0));
   p[1] = std::nextafter(domain.yMin(), real_t(0));
   p[2] = std::nextafter(domain.zMin(), real_t(0));
   mapPointToPeriodicDomain( isPeriodic, domain, p);
   WALBERLA_CHECK_IDENTICAL( p[0], std::nextafter(domain.xMax(), domain.xMin()) );
   WALBERLA_CHECK_IDENTICAL( p[1], std::nextafter(domain.yMax(), domain.yMin()) );
   WALBERLA_CHECK_IDENTICAL( p[2], std::nextafter(domain.zMax(), domain.zMin()) );

   //maxCorner - epsilon -> maxCorner - epsilon
   p[0] = std::nextafter(domain.xMax(), domain.xMin());
   p[1] = std::nextafter(domain.yMax(), domain.yMin());
   p[2] = std::nextafter(domain.zMax(), domain.zMin());
   mapPointToPeriodicDomain( isPeriodic, domain, p);
   WALBERLA_CHECK_IDENTICAL( p[0], std::nextafter(domain.xMax(), domain.xMin()) );
   WALBERLA_CHECK_IDENTICAL( p[1], std::nextafter(domain.yMax(), domain.yMin()) );
   WALBERLA_CHECK_IDENTICAL( p[2], std::nextafter(domain.zMax(), domain.zMin()) );

   return EXIT_SUCCESS;
}

} //namespace testing
} //namespace walberla

int main( int argc, char** argv )
{
   return walberla::testing::main(argc, argv);
}
