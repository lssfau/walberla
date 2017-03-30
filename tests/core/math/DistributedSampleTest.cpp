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
//! \file DistributedSampleTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/DistributedSample.h"
#include "core/math/Random.h"
#include "core/math/Sample.h"
#include "core/mpi/MPIManager.h"

using namespace walberla;

int main(int argc, char * argv[])
{
   debug::enterTestMode();

   auto mpiManager = MPIManager::instance();
   mpiManager->initializeMPI( &argc, &argv );

   WALBERLA_MPI_SECTION() { mpiManager->useWorldComm(); }

   math::Sample sample;
   math::DistributedSample disSample;

   for( uint_t i = 0; i < uint_t(1000); ++i )
   {
      real_t value = math::realRandom( real_c( -100 ), real_c( 100 ) );
      sample.insert( value );
      disSample.insert( value );
   }

   sample.mpiGatherRoot();
   disSample.mpiGatherRoot();

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_CHECK_FLOAT_EQUAL( sample.sum(), disSample.sum() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.min(), disSample.min() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.max(), disSample.max() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.range(), disSample.range() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.mean(), disSample.mean() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.variance(), disSample.variance() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.stdDeviation(), disSample.stdDeviation() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.relativeStdDeviation(), disSample.relativeStdDeviation() );
      WALBERLA_CHECK_EQUAL( sample.size(), disSample.size() );
   }

   sample.clear();
   disSample.clear();

   for( uint_t i = 0; i < uint_t(1000); ++i )
   {
      real_t value = math::realRandom( real_c( 10 ), real_c( 100 ) );
      sample.insert( value );
      disSample.insert( value );
   }

   sample.mpiAllGather();
   disSample.mpiAllGather();

   WALBERLA_CHECK_FLOAT_EQUAL( sample.sum(), disSample.sum() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.min(), disSample.min() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.max(), disSample.max() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.range(), disSample.range() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.mean(), disSample.mean() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.variance(), disSample.variance() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.stdDeviation(), disSample.stdDeviation() );
   WALBERLA_CHECK_FLOAT_EQUAL( sample.relativeStdDeviation(), disSample.relativeStdDeviation() );
   WALBERLA_CHECK_EQUAL( sample.size(), disSample.size() );
   
   sample.clear();
   disSample.clear();
   
   sample.mpiGatherRoot();
   disSample.mpiGatherRoot();
   
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_CHECK_EQUAL( sample.size(), uint_t(0) );
      WALBERLA_CHECK_EQUAL( disSample.size(), uint_t(0) );
   }
   
   sample.clear();
   disSample.clear();   
   
   sample.mpiAllGather();
   disSample.mpiAllGather();
   
   WALBERLA_CHECK_EQUAL( sample.size(), uint_t(0) );
   WALBERLA_CHECK_EQUAL( disSample.size(), uint_t(0) );
   
   if(  mpiManager->numProcesses() > 1 )
   {
      sample.clear();
      disSample.clear();
      
      WALBERLA_NON_ROOT_SECTION()
      {
         for( uint_t i = 0; i < uint_t(1000); ++i )
         {
            real_t value = math::realRandom( real_c( -100 ), real_c( -10 ) );
            sample.insert( value );
            disSample.insert( value );
         }
      }
      
      sample.mpiGatherRoot();
      disSample.mpiGatherRoot();

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_FLOAT_EQUAL( sample.sum(), disSample.sum() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.min(), disSample.min() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.max(), disSample.max() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.range(), disSample.range() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.mean(), disSample.mean() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.variance(), disSample.variance() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.stdDeviation(), disSample.stdDeviation() );
         WALBERLA_CHECK_FLOAT_EQUAL( sample.relativeStdDeviation(), disSample.relativeStdDeviation() );
         WALBERLA_CHECK_EQUAL( sample.size(), disSample.size() );
      }
      
      sample.clear();
      disSample.clear();
      
      WALBERLA_NON_ROOT_SECTION()
      {
         for( uint_t i = 0; i < uint_t(1000); ++i )
         {
            real_t value = math::realRandom( real_c( -100 ), real_c( 100 ) );
            sample.insert( value );
            disSample.insert( value );
         }
      }
      
      sample.mpiAllGather();
      disSample.mpiAllGather();

      WALBERLA_CHECK_FLOAT_EQUAL( sample.sum(), disSample.sum() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.min(), disSample.min() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.max(), disSample.max() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.range(), disSample.range() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.mean(), disSample.mean() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.variance(), disSample.variance() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.stdDeviation(), disSample.stdDeviation() );
      WALBERLA_CHECK_FLOAT_EQUAL( sample.relativeStdDeviation(), disSample.relativeStdDeviation() );
      WALBERLA_CHECK_EQUAL( sample.size(), disSample.size() );
   }
}
