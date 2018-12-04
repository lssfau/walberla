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
//! \file IntersectionRatioTest.cpp
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/geometry/IntersectionRatio.h"
#include "geometry/bodies/Sphere.h"
#include "geometry/bodies/AABBBody.h"

#include <random>


namespace walberla {

void testPointUnitSphere( const math::Vector3<real_t> & p )
{
   static const geometry::Sphere UNIT_SPHERE( math::Vector3<real_t>( 0, 0, 0 ), real_t( 1 ) );
   static const real_t EPSILON = real_t(1e-4);

   real_t q = lbm::intersectionRatioBisection( UNIT_SPHERE, p, -p, EPSILON );
   
   Vector3<real_t> intersectionPoint = p + (-p) * q;
   real_t intersectionRadius = intersectionPoint.length();

   WALBERLA_CHECK_LESS( std::fabs( intersectionRadius - real_t( 1 ) ), EPSILON );

   q = lbm::intersectionRatioSphere( UNIT_SPHERE, p, -p );
   
   intersectionPoint = p + ( -p ) * q;
   intersectionRadius = intersectionPoint.length();
   
   WALBERLA_CHECK_LESS( std::fabs( intersectionRadius - real_t( 1 ) ), EPSILON );
   
   q = lbm::intersectionRatio( UNIT_SPHERE, p, -p, EPSILON );
   
   intersectionPoint = p + ( -p ) * q;
   intersectionRadius = intersectionPoint.length();
   
   WALBERLA_CHECK_LESS( std::fabs( intersectionRadius - real_t( 1 ) ), EPSILON );
}

void testUnitSphere()
{
   testPointUnitSphere( Vector3<real_t>(  1,  1,  1 ) );
   testPointUnitSphere( Vector3<real_t>(  1,  1, -1 ) );
   testPointUnitSphere( Vector3<real_t>(  1, -1,  1 ) );
   testPointUnitSphere( Vector3<real_t>(  1, -1, -1 ) );
   testPointUnitSphere( Vector3<real_t>( -1,  1,  1 ) );
   testPointUnitSphere( Vector3<real_t>( -1,  1, -1 ) );
   testPointUnitSphere( Vector3<real_t>( -1, -1,  1 ) );
   testPointUnitSphere( Vector3<real_t>( -1, -1, -1 ) );
}

void testAABB()
{

   static const math::Vector3<real_t> ZERO( real_t( 0 ), real_t( 0 ), real_t( 0 ) );
   static const math::Vector3<real_t> UNIT( real_t( 1 ), real_t( 1 ), real_t( 1 ) );
   static const real_t EPSILON = real_t(1e-4);

   std::mt19937 randomEngine;

   std::vector<math::AABB> testAABBs;
   testAABBs.emplace_back( -UNIT, UNIT );
   testAABBs.emplace_back(  ZERO, UNIT );
   testAABBs.emplace_back( -UNIT, ZERO );

   for( auto aabbIt = testAABBs.begin(); aabbIt != testAABBs.end(); ++aabbIt )
   {
      const math::AABB outerAABB = aabbIt->getScaled( real_t( 2 ) );
      std::vector< std::pair< Vector3<real_t>, Vector3<real_t> > > testPoints;

      for( int i = 0; i < 100; ++i )
      {
         Vector3<real_t> outerPoint, innerPoint;
         do { outerPoint = outerAABB.randomPoint( randomEngine ); } while( aabbIt->contains( outerPoint ) );
         innerPoint = aabbIt->randomPoint( randomEngine );
         testPoints.emplace_back( outerPoint, innerPoint - outerPoint );
      }
      
      for( auto pointIt = testPoints.begin(); pointIt != testPoints.end(); ++pointIt )
      {
         const Vector3<real_t> & fluidPoint = pointIt->first;
         const Vector3<real_t> & direction  = pointIt->second;

         real_t q = lbm::intersectionRatio( *aabbIt, fluidPoint, direction, EPSILON );
         Vector3<real_t> intersectionPoint = fluidPoint + direction * q;
         WALBERLA_CHECK_LESS( std::fabs( aabbIt->sqSignedDistance( intersectionPoint ) ), EPSILON * EPSILON );

         q = lbm::intersectionRatioBisection( *aabbIt, fluidPoint, direction, EPSILON );
         intersectionPoint = fluidPoint + direction * q;
         WALBERLA_CHECK_LESS( std::fabs( aabbIt->sqSignedDistance( intersectionPoint ) ), EPSILON * EPSILON );
      }
   }

}

int main( int /*argc*/, char ** /*argv*/ )
{   
   testUnitSphere();
   testAABB();

   return EXIT_SUCCESS;
}

} // namespace walberla



int main( int argc, char ** argv ) {

   walberla::debug::enterTestMode();

   return walberla::main( argc, argv );
}



