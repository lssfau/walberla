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
//! \file SampleTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Sample.h"
#include "core/mpi/MPIManager.h"

#include <algorithm>
#include <cmath>
#include <string>


using namespace walberla;
using walberla::math::Sample;

void testOneElement(real_t r)
{
   Sample sample;
   sample.insert(r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.sum(),    r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.min(),    r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.max(),    r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.range(),  real_c(0));

   WALBERLA_CHECK_FLOAT_EQUAL(sample.mean(),   r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.median(), r);

   WALBERLA_CHECK_FLOAT_EQUAL(sample.variance(),     real_c(0));
   WALBERLA_CHECK_FLOAT_EQUAL(sample.stdDeviation(), real_c(0));
   if( std::fabs(r) > 1e-12 )
      WALBERLA_CHECK_FLOAT_EQUAL(sample.relativeStdDeviation(), real_c(0));

   WALBERLA_CHECK_FLOAT_EQUAL( sample.mad(), real_c(0) );

   WALBERLA_CHECK_FLOAT_EQUAL(sample.cummulativeDistributionFunction(r - 1), real_c(0));
   WALBERLA_CHECK_FLOAT_EQUAL(sample.cummulativeDistributionFunction(r + 1), real_c(1));

   WALBERLA_CHECK_FLOAT_EQUAL(sample.quantile(real_c(0)), r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.quantile(real_c(0.5)), r);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.quantile(real_c(1)), r);

}

void testTwoElements( const real_t r1, const real_t r2 )
{
   const real_t minR = std::min( r1, r2 );
   const real_t maxR = std::max( r1, r2 );

   Sample sample;
   sample.insert(r1);
   sample.insert(r2);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.sum(),    r1+r2);
   WALBERLA_CHECK_FLOAT_EQUAL(sample.min(),    minR );
   WALBERLA_CHECK_FLOAT_EQUAL(sample.max(),    maxR );
   WALBERLA_CHECK_FLOAT_EQUAL(sample.range(),  std::fabs(r1-r2));
   WALBERLA_CHECK_FLOAT_EQUAL(sample.mean(),   (r1+r2)/real_c(2));
   WALBERLA_CHECK_FLOAT_EQUAL(sample.median(), (r1+r2)/real_c(2));
   real_t mu = sample.mean();
   WALBERLA_CHECK_FLOAT_EQUAL(sample.variance(), ((r1-mu)*(r1-mu) + (r2-mu)*(r2-mu)) / real_c(2) );
   WALBERLA_CHECK_FLOAT_EQUAL(sample.stdDeviation(), std::fabs(r1 - r2) / real_t(2) );
   if( sample.mean() > 1e-12 )
      WALBERLA_CHECK_FLOAT_EQUAL(sample.relativeStdDeviation(), sample.stdDeviation() / sample.mean() );

   WALBERLA_CHECK_FLOAT_EQUAL( sample.mad(), real_c( std::fabs( r1 - r2 ) / real_c(2) ) );

   if( minR + maxR > real_t(0) )
      WALBERLA_CHECK_FLOAT_EQUAL( sample.giniCoefficient(), real_t(2) * maxR / (minR + maxR) - real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL(sample.cummulativeDistributionFunction( minR - real_c(1)), real_c(0));
   if(r1 > r2 || r1 < r2)
      WALBERLA_CHECK_FLOAT_EQUAL(sample.cummulativeDistributionFunction(real_c(0.5) * (r2+r1)), real_c(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(sample.cummulativeDistributionFunction( maxR + real_c(1)), real_c(1));

   WALBERLA_CHECK_FLOAT_EQUAL(sample.quantile(real_c(0)), minR );
   WALBERLA_CHECK_FLOAT_EQUAL(sample.quantile(real_c(1)), maxR );
}

void testStaticSample()
{
   Sample sample;
   sample.insert( real_c(-17) );
   sample.insert( real_c(0) );
   sample.insert( real_c(2) );
   sample.insert( real_c(9) );
   sample.insert( real_c(12) );
   sample.insert( real_c(21) );
   sample.insert( real_c(36) );
   sample.insert( real_c(100) );
   sample.insert( real_c(112) );

   WALBERLA_CHECK_FLOAT_EQUAL( sample.min(), real_c(-17) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.max(), real_c(112) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.range(), real_c(112 - -17) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.median(), real_c(12) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.mean(), real_c(30.5555555555556000) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.variance(), real_c(1821.8024691358000000) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.stdDeviation(), real_c(42.6825780516571000) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.mad(), real_c(12) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.quantile( real_c(0.25) ), real_c(2) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.quantile( real_c(0.5) ),  real_c(12) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.quantile( real_c(0.75) ), real_c(36) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.cummulativeDistributionFunction( real_c(2) ),  real_c( 0.33333333333333333333333333333333 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.cummulativeDistributionFunction( real_c(12) ), real_c( 0.55555555555555555555555555555556 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.cummulativeDistributionFunction( real_c(36) ), real_c( 0.77777777777777777777777777777778 ) );

   sample.insert( real_c( 150 ) );

   WALBERLA_CHECK_FLOAT_EQUAL( sample.median(), real_c(16.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.mad(), real_c(18) );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.giniCoefficient(), real_c( 0.749542483660131 ) );

   std::vector<int> intValues;
   intValues.push_back( -17 );
   intValues.push_back( 0 );
   intValues.push_back( 2 );
   intValues.push_back( 9 );
   intValues.push_back( 12 );
   intValues.push_back( 21 );
   intValues.push_back( 36 );
   intValues.push_back( 100 );
   intValues.push_back( 112 );
   intValues.push_back( 150 );

   Sample sample2;

   sample2.castToRealAndInsert( intValues.front() );
   sample2.castToRealAndInsert( sample2.end(), intValues[1] );
   sample2.castToRealAndInsert( intValues.begin() + 2u, intValues.end() );

   WALBERLA_CHECK_EQUAL( sample, sample2 );
}

Sample makeUniformDistributedSample(size_t n)
{
   Sample statReal;

   for(size_t i = 1; i <= n; ++i)
   {
      statReal.insert(real_c(i));
   }

   return statReal;
}

void testUniformDistributedSample(const Sample & uniformDist)
{
   size_t n = uniformDist.size();
   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.sum(),    real_c(n*(n+1)/2));
   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.min(),    real_c(1));
   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.max(),    real_c(n));

   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.mean(),   real_c(n+1) / real_c(2));
   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.median(), ( n % 2 ) ? real_c(n/2+1) : real_c(n/2 + n/2 + 1) / real_c(2));

   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.variance(),     real_c(n*n-1)/real_c(12) );
   WALBERLA_CHECK_FLOAT_EQUAL(uniformDist.stdDeviation(), std::sqrt(uniformDist.variance()));
}

void testNoMPI()
{
   for(int i = -5; i <= 5; ++i)
   {
      testOneElement(real_c(i));
   }

   for(int i1 = -5; i1 <= 5; ++i1)
      for(int i2 = -5; i2 <= 5; ++i2)
      {
         testTwoElements(real_c(i1), real_c(i2));
      }

   testStaticSample();
}

void testAllGather(Sample processDist)
{
   size_t numProcesses = numeric_cast<size_t>(MPIManager::instance()->numProcesses());

   size_t sizeOneProcess = processDist.size();

   processDist.mpiAllGather();

   WALBERLA_CHECK_EQUAL(processDist.size(), sizeOneProcess * numProcesses);
   testUniformDistributedSample(processDist);
}

void testGather(Sample processDist)
{
   int numProcesses = MPIManager::instance()->numProcesses();

   size_t sizeOneProcess = processDist.size();

   for(int gatherRank = 0; gatherRank < numProcesses; ++gatherRank)
   {
      processDist.mpiGather(gatherRank);
      WALBERLA_EXCLUSIVE_SECTION(gatherRank)
      {
         WALBERLA_CHECK_EQUAL(processDist.size(), sizeOneProcess * uint_c(numProcesses));
         testUniformDistributedSample(processDist);
      }
      else
      {
         WALBERLA_CHECK(processDist.empty());
      }
   }
}

void testGatherRoot(Sample processDist)
{
   size_t numProcesses = numeric_cast<size_t>(MPIManager::instance()->numProcesses());

   size_t sizeOneProcess = processDist.size();

   processDist.mpiGatherRoot();
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_CHECK_EQUAL(processDist.size(), sizeOneProcess * numProcesses);
      testUniformDistributedSample(processDist);
   }
   WALBERLA_NON_ROOT_SECTION()
   {
      WALBERLA_CHECK(processDist.empty());
   }
}

void testMPI()
{
   size_t numProcesses = numeric_cast<size_t>( MPIManager::instance()->numProcesses() );
   size_t rank         = numeric_cast<size_t>( MPIManager::instance()->rank() );

   for( size_t n = numProcesses; n < numProcesses * 20; ++n )
   {
      Sample dist = makeUniformDistributedSample(n);
      testUniformDistributedSample(dist);

      Sample processDist;
      auto first = dist.begin();
      std::advance(first, numeric_cast< Sample::difference_type >(rank * (n / numProcesses)));
      auto last = first;
      std::advance(last, numeric_cast< Sample::difference_type >((n / numProcesses)));

      for( ; first != last; ++first )
         processDist.insert(*first);

      WALBERLA_CHECK_EQUAL(processDist.size(), n / numProcesses);

      testAllGather(processDist);
      testGather(processDist);
      testGatherRoot(processDist);
   }
}

int main(int argc, char * argv[])
{
   debug::enterTestMode();

   auto mpiManager = MPIManager::instance();
   mpiManager->initializeMPI( &argc, &argv );

   WALBERLA_MPI_SECTION() { mpiManager->useWorldComm(); }

   WALBERLA_ROOT_SECTION()
   {
      testNoMPI();
   }

   testMPI();
}
